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

Class
    Orthog

Description
    Modified Chebyshev algorithm algorithm for moment inversion

    @article{Gautschi1994,
        author = {Gautschi, Walter},
        title = {Algorithm 726: \mbox{ORTHPOL} --- \mbox{A} Package of Routines for
                generating Orthogonal Polynomials and Gauss-Type Quadrature Rules},
        journal = {ACM Trans. Math. Softw.},
        volume = {20},
        number = {1},
        year = {1994}
    }

SourceFiles
    Orthog.C

Author
    LTFD. All rights reserved

    Laboratory of Thermo-Fluid Dynamics (LTFD)
    Chemical Engineering Program (PEQ)
    Alberto Luiz Coimbra Institute for Graduate Studies and Research in Engineering (COPPE)
    Federal University of Rio de Janeiro Rio de Janeiro, Brazil

    Code implementation:
    Jovani L. Favero

\*---------------------------------------------------------------------------*/

#ifndef Orthog_H
#define Orthog_H

#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "EigenSolver.H"
#include "volFieldsFwd.H"
#include "RectangularMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Orthog Declaration
\*---------------------------------------------------------------------------*/

class Orthog
{
    // Private data

        //- Weights
        scalarField weights_;

        //- Abscissas
        scalarField abscissas_;

        //- Error code return
        label Error_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Orthog(const Orthog&);

        //- Disallow default bitwise assignment
        void operator=(const Orthog&);


public:

    //- Runtime type information
    TypeName("Orthog");


    // Constructors

        //- Construct from scalar field moments
        Orthog(const scalarField& nMom);


    // Destructor

        virtual ~Orthog()
        {}

    // Member Functions

        //- Calculate factorial
        scalar factorial(label n);

        //- Return weights
        const scalarField& weights() const
        {
            return weights_;
        }

        //- Return abscissas
        const scalarField& abscissas() const
        {
            return abscissas_;
        }

        //- Return Error code
        const label& Error() const
        {
            return Error_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
