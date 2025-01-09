# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Mon 30 Dec 2024 02:28:48



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
               value = 125.11,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

Me = Parameter(name = 'Me',
               nature = 'external',
               type = 'real',
               value = 0.000511,
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
               value = 1,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '0.30282212',
               texname = '\\text{ee}')

lambdah = Parameter(name = 'lambdah',
                    nature = 'internal',
                    type = 'real',
                    value = '0.0057',
                    texname = '\\text{lambdah}')

lambdacc = Parameter(name = 'lambdacc',
                     nature = 'internal',
                     type = 'real',
                     value = '5',
                     texname = '\\text{lambdacc}')

