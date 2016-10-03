/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include <cmath>
#include <cfloat>
#include "Orthog.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    extern "C"
    {
        void dcheb_
                (
                    label *, scalar *, scalar *, scalar *, scalar *, scalar *,
                    scalar *, label *, scalar *, scalar *, scalar *, scalar *, scalar *
                );
        void dgauss_
                (
                    label *, scalar *, scalar *, scalar *, scalar *, scalar *,
                    label *, scalar *
                );

    }

    defineTypeNameAndDebug(Orthog, 0);

    scalar Orthog::factorial(label n)
    {
        return n >= 1 ? n * factorial(n-1) : 1;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given a scalar field
Foam::Orthog::Orthog
(
    const scalarField& moments
)
:
    weights_(moments.size()/2),
    abscissas_(moments.size()/2)
{
    //- Number of moments to deal
    label chebN, chebErr, gaussErr, nMom, chebNo;

    scalar
        *chebA, *chebB, *chebMu, *chebAlpha, *chebBeta,
        *chebS, *chebS0, *chebS1, *chebS2, *absc, *weight, *e,
         eps, minf, maxf;

    nMom = moments.size();
    chebN = nMom/2;
    chebNo = chebN;
    chebErr = -1;
    gaussErr = -1;

    //- Define machine precision, originaly
    // done by d1mach_ subroutine from orthog package
    eps = DBL_EPSILON/FLT_RADIX;
    minf = DBL_MIN;
    maxf = DBL_MAX;

    //- Allocate memory and test the allocation
    chebA     = new scalar[nMom-1];
    chebB     = new scalar[nMom-1];
    chebMu    = new scalar[nMom];
    chebAlpha = new scalar[chebN];
    chebBeta  = new scalar[chebN];
    chebS     = new scalar[chebN];
    chebS0    = new scalar[nMom];
    chebS1    = new scalar[nMom];
    chebS2    = new scalar[nMom];
    absc      = new scalar[chebN];
    weight    = new scalar[chebN];
    e         = new scalar[chebN];

    if (
            !chebA || !chebB || !chebMu || !chebAlpha || !chebBeta || !chebS
        ||
            !chebS0 || !chebS1 || !chebS2 || !absc || !weight || !e
       )
    {
        FatalErrorIn("Orthog::Orthog::(..)")
            << "\nError: not enough memory for allocation!"
            << abort(FatalError);
    }

    //- For standard moments. It can be improved to generalized moments handling
    for(label i=0; i<nMom-1; i++)
    {
        chebA[i] = 0.0;
        chebB[i] = 0.0;
    }

    forAll(moments, i)
    {
        chebMu[i] = moments[i];
    }

    //- If mu0 is small, calculate gamma distribution  standard moments with zero weight
//     if(chebMu[0] < 100*minf)
//         for(label i=0; i<nMom; i++)
//             chebMu[i] = factorial(i)*SMALL;


    Error_ = 0;

    //- Initializing weight and abscissas
    for(label i=0; i<chebNo; i++)
    {
        weight[i] = 0.0;
        absc[i]   = 0.0;
    }

//      if( (min(moments)/(max(moments)+minf)) > 0.9999999999  )   Error_ = 1;
    //- Decrease number of quadrature points untill the moment set is realizable
    do
    {

        if(chebN != chebNo)
        {
            Info << "Non realizable moment set " << moments
            << ". Decreasing number of quadrature points to "
            << chebN << endl;
        }

        //- Initializing nulls alphas and betas
        for(label i=0; i<chebN; i++)
        {
            chebAlpha[i] = 0.0;
            chebBeta[i]  = 0.0;
        }


        if(chebMu[0] > SMALL) //100*minf)
        {
            //- Fortran routine to calculate alpha and beta coefficients
            //(chebAlpha, chebBeta) to be used in the dgauss routine
            dcheb_
            (
                &chebN, chebA, chebB, chebMu, chebAlpha, chebBeta, chebS,
                &chebErr, chebS0, chebS1, chebS2, &minf, &maxf
            );

            //- Test for the chebErr return
            switch (chebErr)
            {
                case 0:
                    //Info << "\nchebErr = " << chebErr
                    //        << ". dcheb: normal return." << endl;
                break;
                case 1:
                    Info << "\nchebErr = " << chebErr
                        << ". dcheb: Moment value is lower"
                        <<" than the machine zero" << endl;
                break;
                case 2:
                    Info << "\nchebErr = " << chebErr
                        << ". dcheb: n is not in range." << endl;
                break;
                default:
                    FatalErrorIn("Orthog::Orthog::dcheb_(..)") << "\nchebErr = "
                        << chebErr << ". dcheb: if + is about to underflow if - is"
                        <<" about to overflow."
                        << abort(FatalError);
            }

            //- Fortran routine to calculate weight and abscissas
            dgauss_
            (
                &chebN, chebAlpha, chebBeta, &eps, absc, weight,
                &gaussErr,e
            );

            //- Test for the gaussErr return
            switch (gaussErr)
            {
                case 0:
                    Error_ = 0;
                    //Info << "\ngaussErr = " << gaussErr
                    //        << ". dgauss: normal return." << endl;
                break;
                case -1:
                    Info << "\ngaussErr = " << gaussErr
                        << ". dgauss: chebN is not in range." << endl;
                    Error_ = 1;
                break;
                case -2:
                    Info << "\ngaussErr = " << gaussErr
                        << ". dgauss: one of the beta's"
                        <<" coefficients is negative." << endl;
                    Error_ = 1;
                break;
                default:
                    FatalErrorIn("Orthog::Orthog::dgauss_(..)")
                        << "\ngaussErr = "
                        << gaussErr
                        << ". dgauss: the QR algorithm does not converge within"
                        <<" 30 iterations on evaluating the"
                        << gaussErr << " th eigenvalue."
                        << abort(FatalError);
                    Error_ = 1;
            }

        }

        for(label i=0; i<chebN; i++)
            if(absc[i] < 0.0)
                absc[i] = 0.0;

        chebN--;

    }while((Error_ != 0) && (chebN > 0));

    //- Returning weights and abscissas in decreasing order
    for(label i=0; i<chebNo; i++)
    {
         weights_[i]   =  weight[i];//weight[chebNo-1-i];
         abscissas_[i] =  absc[i];//absc[chebNo-1-i];
    }

    //- Free allocated memory
    delete [] chebA;
    delete [] chebB;
    delete [] chebMu;
    delete [] chebAlpha;
    delete [] chebBeta;
    delete [] chebS;
    delete [] chebS0;
    delete [] chebS1;
    delete [] chebS2;

    delete [] absc;
    delete [] weight;
    delete [] e;

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// ************************************************************************* //
