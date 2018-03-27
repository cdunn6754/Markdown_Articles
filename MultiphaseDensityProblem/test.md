# Multiphase soot density treatment in OpenFOAM 5.x


## Overview
With the twoEquationSoot (tes) soot model the soot is essentially treated as a gas phase specie and we use the graphite janaf data to represent it. This makes it much easier to use the devolatilization and (most of the) thermodynamic/transport machinery that is built in to OF. One significant detriment to this approach is in the calculation of the gas phase density.

The density of the gas phase mixture (at least with the thermo classes used in coalChemistryFoam and its derivatives) is calculated as

$$ \rho = \psi * P $$
in the `psiThermo::rho()` function.

Above, $\psi$ is the fluid compressibility field and $P$ is the pressure field. The pressure field is calculated with the Poisson equation but $\psi$ is determined from a mass weighted average of the individual specie compressibilities ($\psi_i$). This is inappropriate for our soot specie because physically it is small discrete particles of soot and not individual carbon atoms floating in the gas mixture as the use of the graphite janaf chemistry data implies.

The proper multiphase treatment is to calculate the gas density as

$$ \rho = (1 - V_s) \cdot (\psi_{*} * P) + (V_s) \cdot (\rho_{soot})
  \tag{1}$$

Where $V_s$ is the soot volume fraction, $\psi_{*}$ is the fluid compressibility calculated as before but neglecting the contribution of the soot specie. The density of coal soot is empirically known to be
$2000 \left[kg/m^3\right]$.

## Relevant classes and functions in OpenFOAM

The thermo class used in coalChemistryFoam, i.e.  what the `thermo` pointer is pointed at is

```c++
hePsiThermo<psiReactionThermo,
            SpecieMixture<reactingMixture<gasHthermoPhysics>>
            >
```
and gasHthermoPhysics is a typedef for another class

```c++
sutherlandTransport<species::thermo<janafThermo<perfectgas<specie>>,
                                    sensibleEnthalpy
                                   >
                   >
```

Obviously that really becomes confusing quickly but there are only a few functions that are relevant to what we are trying to accomplish.

- `psiThermo::rho()`

  Same function as mentioned above. It calculates the density `volScalarField` based on the multiplication of stored values of inherited member variables `psi_` and `p_`.

- `multiComponentMixture::cellMixture(celli)`

  `multiComponentMixture` is inherited by the `reactingMixture` class. The function calculates the thermodynamic properties for a cell based on the current mass fractions of species in that cell.

  ```c++
      mixture_ = Y_[0][celli]*speciesData_[0];

      for (label n=1; n<Y_.size(); n++)
      {
        mixture_ += Y_[n][celli]*speciesData_[n];
      }
  ```

  Here `speciesData_` is a list of `gasHThermoPhysics` objects, one object entry per specie. For example to find the molecular wieght of $CH_4$, if it were the second specie in the list one could use `speciesData_.W(1)`. The overloaded `+=` operator is used in coordination with the species mass fractions, `Y_`, to calculate mixture quantities. Importantly to us, the compressibility $\psi$, is one of these averaged quantities. We will therefore need to modify this function to exclude soot when taking the mass wieghted average.

- `multiComponentMixture::patchFaceMixture(patchi, facei)`

  Same as `multiComponentMixture::cellMixture(celli)` but determines the mixture thermo properties on a boundary face, not at a cell center. Everything that we do later to modify `cellMixture` needs to be done in this function as well.

- `hePsiThermo::calculate()`

  This function loops over all cells and calls `multiComponentMixture::cellMixture(celli)` on each of them. Once the mixture properties are calculated for a cell it sets the main thermodynamic fields of that cell based on the mixture, these include the temperature, compressibility, dynamic viscosity and thermal diffusivity.

  It subsequently calls `multiComponentMixture::patchFaceMixture(celli)` for each boundary face to similarly calculate and set the thermodynamic properties on the boundaries.

  It is in this function that we will need to differentiate between the standard `cellMixture` and a `sootCellMixture` function, which excludes soot from the mixture, when setting the values for the `psi_` field.

## Source code modifications
It turns out that all of the needed modifications can be performed at the level of `hePsiThermo` which is particularly convenient since that is one of the classes that can be specified at run-time in the thermophysicalProperties dictionary file. Assuming the appropriate modifications can be made and compiled in a library with a renamed class, here it is called `sootHePsiThermo`, we can just select that class in the dictionary and avoid further changes to the solver source code. Unfortunately as described below there is an additional complication that is resolved with a solver modification but at least major changes to the compile time solver models are avoided.

The overall plan is as follows:
1. Copy hePsiThermo.H and hePsiThermo.C to a user source directory and rename the files and class as `sootHePsiThermo`
2. Write, within the new class, a `sootCellMixture(celli)` function that excludes soot from consideration
3. Modify the `calculate()` function to use the new `sootCellMixture(celli)` function when setting global field `psi_`
4. Add `sootVolume_` field and `updateSootVolume()` function to class
5. Create a function `sootHePsiThermo::rho()` that overrides `psiThermo::rho()` and uses equation 1 to calculate the density
6. Modify solver to downcast the `thermo` pointer, allowing us direct access to the functions we just wrote


#### 1. Copy and rename the class
Copy the class from the `thermophysicalModels/basic/psiThermo/` directory by copying hePsiThermo.H and hePsiThermo.C to a user source directory. You will also need to copy the source file in which the new type will be created in the runtime selection table. That file is named psiReactionThermos.C and is located at `thermophysicalModels/reactionThermo/psiReactionThermo/`. Then change the file names (I used 'sootHePsiThermo') and use sed to change all 'hePsiThermo' to 'sootHePsiThermo' e.g.

```
sed -i -e 's/hePsiThermo/sootHePsiThermo/g' *
```

Be sure not to change any other names, it is easy to get carried away and accidentally change 'psiReactionThermo' to 'sootPsiReactionThermo'for instance.

You will also need to copy a Make directory if you want to use wmake, which I would recommend. I took it from the `thermophysicalModels/reactionThermo` directory. The only file that you will be compiling (i.e. listing in the `Make/files` file) is sootPsiThermos.C, the `Make/options` file should be okay as is but if when compiling you get an error just be sure to add whichever library/header file it complains about missing to it.

#### 2. Writing sootCellMixture function
As mentioned above the built in version of this function, `multicomponentMixture::cellMixture(celli)`, needs to be modified to exclude the soot specie from consideration. Here is the implementation of the new function

```c++
template<class BasicPsiThermo, class MixtureType>
typename MixtureType::thermoType
Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::sootCellMixture
(
    const label celli
) const
{

    // Get those mass fractions
    const PtrList<volScalarField>& Y_ = MixtureType::Y();

    // Hope that SOOT isn't the first specie listed
    typename MixtureType::thermoType mixture =
        Y_[0][celli]*MixtureType::speciesData()[0];

    // The whole point is to avoid SOOT
    if (Y_[0].name() != "SOOT")
    {
        for (label n=1; n<Y_.size(); n++)
        {
            if (Y_[n].name() != "SOOT")
            {
                mixture += Y_[n][celli]*MixtureType::speciesData()[n];
            }
        }
    }
    else
    {
        // If SOOT was the first one then reset and loop through the rest
        mixture = Y_[1][celli]*MixtureType::speciesData()[1];

        for (label n=1; n<Y_.size(); n++)
        {
            mixture += Y_[n][celli]*MixtureType::speciesData()[n];
        }
    }

    return mixture;
}
```

Fortunately we have the `MixtureType` template arguement here that allows us to refer to anything we need from the mixture classes so this can almost be directly copied from the main OpenFOAM implementation in the `multiComponentMixture` class.

Without going into the details here, you will need to create another function (mine is named `sootPatchFaceMixture(patchi,facei)`) that does for the boundary faces what `sootCellMixture(celli)` does for the cell centers (exclude soot from the calculation). The necessary code changes are entirely analogous.

#### 3. Modify `calculate()` to use our new `sootCellMixture(celli)`
Now that we can calculate important thermodynamic variables in a cell while excluding soot we need to modify this function to set the main themo member `volScalarField psi_` accordingly. Here is the modified implementation which utilized `sootCellMixture(celli)` to calculate $\psi$.

```c++
template<class BasicPsiThermo, class MixtureType>
void Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        // Use the sootCellMixture function to exclude
        // soot from the calculation
        const typename MixtureType::thermoType sootMixture_ =
            this->sootCellMixture(celli);

        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        psiCells[celli] = sootMixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& sootMixture_ =
                    this->sootPatchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = sootMixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& sootMixture_ =
                    this->sootPatchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                ppsi[facei] = sootMixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}
```

#### 4. Add a `sootVolume_` field and `updateSootVolume()` function to the class

In order to calculate the density in accordance with equation 1 we need to determine the soot volume fraction field. Since this field is of interest in its own right another valid approach might be to create the field within the solver and then just pass it to this class to make the density calculation. Here I have made it a member of the class and just added an acess function for it.

You can see the addition of the field on my [github](https://github.com/cdunn6754/OpenFOAM_5.x_Libraries), there is nothing too special about it. For convenience I also added a soot density member variable, `sootDensity_`, that is hardcoded to $2000 \, [kg/m^3]$.

The function `updateSootVolume()` is taken with only minor modifications from the 'greyMeanSolidAbsoprtionEmission' radiation absorption/emission model. I thought it was a little complicated when I first saw it but I now think it is the best way to calculate volume fraction of a specie. The basic idea is to first calculate something like the mixture specific volume

$$
  \nu_{mix}
  =
  \sum_{i}^{species} Y_i / \rho_i(T,P)
$$

Where $Y_i$ is a specie mass fraction and the specie density $\rho_i$ is calculated with the ideal gas law for the individual specie. We sum $Y_i / \rho_i$ which has units $[V_i / \text{unit mass}\,]$. Summing these yields the total volume of all gas species per unit mass ($\dot{=} \, V_{mix} / \text{unit mass}$). Then taking $\nu_{soot} = Y_{soot} / \rho_{soot}$ we can calculate the volume fraction as $\nu_{soot} / \nu_{mix}$ which has units of $V_{soot} / V_{mix}$. It is important here to use the known density of soot (2000 $[kg/m^3]$) rather than the ideal gas law prediction based on janaf thermodynamic data.

A problem that remains in this approach is the question of what to do with the boundary values. The density field is a `volScalarField` and it therefore needs to have boundary conditions defined and we will in turn need boundary values for the soot volume fraction to calculate it. I think it should be possible to calculate boundary face valus since the `p_`, `T_` and `Y` fields are also `volScalarField`s but for now I am going to assume zero soot volume fraction at the boundaries (the last line in the function). At any rate if you do nothing the uninitialized boundary values will cause drastic density fluctuations and crash the simulation immediatly.


```c++
template<class BasicPsiThermo, class MixtureType>
void Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::updateSootVolume()
{

    // Hardcoded soot density
    scalar sootDensity(2000.0); // kg\m^3 from dasgupta thesis

    // To be the sum overall species of [m^3_specie / kg_total]
    scalarField specificVolumeSum =
    scalarField(this->sootVolume_.size(), 0.0);

    // As we iterate we will grab the SOOT specie index
    label sootIdx(-1);

    // Pointer to the mixture for this thermo
    basicSpecieMixture& mixture_ = this->composition();

    forAll(this->Y(), specieI)
    {
        const scalarField& Yi = mixture_.Y()[specieI];
        const word specieName = mixture_.Y()[specieI].name();

        if (specieName == "SOOT")
        {
            sootIdx = specieI;
            specificVolumeSum += Yi/sootDensity;
        }
        else
        {
            // loop through cells for non-constant density
            forAll(specificVolumeSum, celli)
            {
                specificVolumeSum[celli] += Yi[celli]/
                    mixture_.rho(specieI, this->p_[celli], this->T_[celli]);
            }
        }

    }// end loop through species

    // now find and set the soot volume fraction as
    // [V_soot/kg_total] / [V_total/kg_total]
    this->sootVolume_.primitiveFieldRef() =
        (this->Y()[sootIdx]/sootDensity) / (specificVolumeSum);

    // Set this to 0 so that psi_* P_ is used for the boundary density field.
    this->sootVolume_.boundaryFieldRef() = 0.0;
}
```

#### 5. Create the `sootHePsiThermo::rho()` function

Finally we are ready to actually write a new density calculation function. All of the work is already done and we can just write out the density function as described in equation 1 using the member variables we added earlier.

```c++
template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::sootHePsiThermo<BasicPsiThermo, MixtureType>::rho() const
{

    return (1.0 - this->sootVolume_)*(this->p_*this->psi_) +
        (this->sootVolume_) * this->sootDensity_;
}
```

#### 6. Modifiy the solver to allow access to new class

As mentioned above I thought this would be a good approach because the use of `sootHePsiThermo` rather than `hePsiThermo` is specified at runtime from the 'thermophysicalProperties' dictionary. That means that no solver modifications are necessary. Unfortunately I was mistaken because the thermo class within the solver is not created directly but instantiated within the combusition class. The combustion class then passes an upcasted reference to the solver, here is the code from the coalChemistryFoam createFields.H file

```c++
Info<< "Creating combustion model\n" << endl;

autoPtr<combustionModels::psiCombustionModel> combustion
(
    combustionModels::psiCombustionModel::New(mesh)
);

psiReactionThermo& thermo = combustion->thermo();
```

Examining the thermo class in a little more detail we can see why this is possible. First here is the class again, as shown at the top of this document I have replace the specific template parameters used with more general names (those used for the template parameters in the source) for brevity.

```c++
hePsiThermo<BasicThermo,MixtureType>
```

You can look at the documentaton of the classes and discover that while `hePsiThermo` does not directly inherit from its template parameters, `BasicThermo` and `MixtureType`, it does inherit from `heThermo<BasicThermo,MixtureType>`. And `heThermo<BasicThermo,MixtureType>` inherits from both `BasicThermo` and `MixtureType`. So in the specific case where we `BasicThermo` is `psiReactionThermo` we know that `hePsiThermo` indirectly inherits from `psiReactionThermo` and it can therefore be upcast as implied in the code snippet above. I think they actually do a `dynamic_cast` on the original thermo reference within the combustion model to upcast and then just pass that member reference variable here.

The problem with the upcast to `psiReactionThermo` is that all of the functions we just wrote in `hePsiThermo` are now inaccessible (you can't use a base class reference to access derived class functions when those functions they aren't present in the base class, and even then only if they are virtual functions). So the options are to either create a new combustion model too, with only the type of the thermo reference changed. Or to just use a downcast within the solver to change the reference to `psiReactionThermo` to a reference to `sootHePsiThermo`. The second option is used here.

Here is my new version of the createFields file for the SootCoalFoam solver in which the thermo pointer is downcast using `dynamic_cast`

```c++
autoPtr<combustionModels::psiCombustionModel> combustion
(
    combustionModels::psiCombustionModel::New(mesh)
);

psiReactionThermo& baseThermo = combustion->thermo();


// Downcast baseThermo to thermo
// changes type from psiReactionThermo to sootHePsiThermo
// which enables us to use the functions in sootHePsiThermo
sootHePsiThermo<
    psiReactionThermo,
    SpecieMixture< reactingMixture< gasHThermoPhysics > >
    > &thermo =
    dynamic_cast<
        sootHePsiThermo<
            psiReactionThermo,
            SpecieMixture<reactingMixture<gasHThermoPhysics > >
            > &>
                                        (baseThermo);
```
It is really messy, maybe some typdefs would be helpful to understand here but all we are doing is telling the pointer that it now points to type

```c++
hePsiThermo<psiReactionThermo,
            SpecieMixture<reactingMixture<gasHthermoPhysics>>
            >
```

rather than just `psiReactionThermo`

I assume that the creators of OpenFOAM wrote it that way for a reason and I'm a little afraid that having the thermo pointer like that will cause a problem but so far in my testing I have not encountered any problems.

You will need to link to the new library we created with the `sootHePsiThermo` class before compiling the solver. You will also need to include some header files for the additional classes needed in the downcast like `SpecieMixture`, `gasHthermoPhysics` and `reactingMixture`.
