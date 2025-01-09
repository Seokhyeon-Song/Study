# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Sat 28 Dec 2024 18:41:08



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# This is a default parameter object representing the renormalization scale (MU_R).
MU_R = Parameter(name = 'MU_R',
                 nature = 'external',
                 type = 'real',
                 value = 91.188,
                 texname = '\\text{\\mu_r}',
                 lhablock = 'LOOP',
                 lhacode = [1])

# User-defined parameters.
Mh = Parameter(name = 'Mh',
               nature = 'external',
               type = 'real',
               value = 125.11,
               texname = '\\text{Mh}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

Mcc = Parameter(name = 'Mcc',
                nature = 'external',
                type = 'real',
                value = 60,
                texname = '\\text{Mcc}',
                lhablock = 'MASS',
                lhacode = [ 9000001 ])

Mele = Parameter(name = 'Mele',
                 nature = 'external',
                 type = 'real',
                 value = 0.000511,
                 texname = '\\text{Mele}',
                 lhablock = 'MASS',
                 lhacode = [ 11 ])

Mmu = Parameter(name = 'Mmu',
                nature = 'external',
                type = 'real',
                value = 0.10566,
                texname = '\\text{Mmu}',
                lhablock = 'MASS',
                lhacode = [ 13 ])

Mtau = Parameter(name = 'Mtau',
                 nature = 'external',
                 type = 'real',
                 value = 1.777,
                 texname = '\\text{Mtau}',
                 lhablock = 'MASS',
                 lhacode = [ 15 ])

Wh = Parameter(name = 'Wh',
               nature = 'external',
               type = 'real',
               value = 1,
               texname = '\\text{Wh}',
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

