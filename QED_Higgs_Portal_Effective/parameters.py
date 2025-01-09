from object_library import all_parameters, Parameter
from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 125.,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

Me = Parameter(name = 'Me',
               nature = 'external',
               type = 'real',
               value = 0.0005110000000000001,
               texname = '\\text{Me}',
               lhablock = 'MASS',
               lhacode = [ 11 ])

MM = Parameter(name = 'MM',
               nature = 'external',
               type = 'real',
               value = 0.10566,
               texname = '\\text{MM}',
               lhablock = 'MASS',
               lhacode = [ 13 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 6.38233934e-03,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 132.50698,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

lambdah = Parameter(name = 'lambdah',
                    nature = 'external',
                    type = 'real',
                    value = '0.00575',
                    texname = '\\lambda',
                    lhablock = 'YUKAWA',
                    lhacode = [ 1 ])

ch = Parameter(name = 'ch',
                    nature = 'external',
                    type = 'real',
                    value = '0.04',
                    texname = '\\text{c}_{h}',
                    lhablock = 'YUKAWA',
                    lhacode = [ 2 ])

mchi = Parameter(name = 'mchi',
                    nature = 'external',
                    type = 'real',
                    value = '62',
                    texname = '\\text{m}_{\\chi}',
                    lhablock = 'YUKAWA',
                    lhacode = [ 3 ])
                    
aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\text{aEW}')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')