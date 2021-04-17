```
dynamicSmagorinsky - Implementation of the dynamic Smagorinsky SGS model
                     as proposed by Lilly (1992) for OpenFOAM.

Copyright Information
    Copyright (C) 1991-2009 OpenCFD Ltd.
    Copyright (C) 2010-2021 Alberto Passalacqua

License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Description
    Implementation of the dynamic Smagorinsky model with coefficients cD and
    cI computed as local average of their face values to avoid numerical
    instabilities.
    Negative values of the effective viscosity are removed by clipping it to
    zero (nuSgs is clipped to -nu)

Target platform
    The code is known to work with OpenFOAM v2012.

Author
    Alberto Passalacqua <apcfd@outlook.com>
    Website: http://www.albertopassalacqua.com

```

# Instructions


1. Download the source code using git:

    git clone git://github.com/AlbertoPa/dynamicSmagorinsky.git

2. Enter the directory where the source code has been extracted, and compile
   it by typing: wmake libso

3. Add the following line to the controlDict of your case:

    libs ( "libOpenFOAM.so" "libdynamicSmagorinskyModel.so" ) ;

4. Specify

    LESModel        dynamicSmagorinsky;
    delta           cubeRootVol;

   in turbulentProperties.

5. Add the subdictionary

    dynamicSmagorinskyCoeffs
    {
      filter    simple;
      ce        1.048;
    }

   to turbulentProperties.

# Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, the producer
of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD®  trade marks.

Detailed information on the OpenFOAM trademark can be found at

 - http://www.openfoam.com/legal/trademark-policy.php
 - http://www.openfoam.com/legal/trademark-guidelines.php

For further information on OpenCFD and OpenFOAM, please refer to

 - http://www.openfoam.com