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

You can see the addition of the field on my [github](https://github.com/cdunn6754/OpenFOAM_5.x_Libraries), there is nothing too special about it. The function `updateSootVolume()` is taken with only minor modifications from the 'greyMeanSolidAbsoprtionEmission' radiation absorption/emission model. I thought it was a little complicated when I first saw it but I now think it is the best way to calculate volume fraction of a specie.
