# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Fri 20 Dec 2024 11:45:40



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
Mh = Parameter(name = 'Mh',
               nature = 'external',
               type = 'real',
               value = 125.11,
               texname = '\\text{Mh}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

ME = Parameter(name = 'ME',
               nature = 'external',
               type = 'real',
               value = 0.000511,
               texname = '\\text{ME}',
               lhablock = 'MASS',
               lhacode = [ 11 ])

MM = Parameter(name = 'MM',
               nature = 'external',
               type = 'real',
               value = 0.10566,
               texname = '\\text{MM}',
               lhablock = 'MASS',
               lhacode = [ 13 ])

ML = Parameter(name = 'ML',
               nature = 'external',
               type = 'real',
               value = 1.777,
               texname = '\\text{ML}',
               lhablock = 'MASS',
               lhacode = [ 15 ])

Mf = Parameter(name = 'Mf',
               nature = 'external',
               type = 'real',
               value = 60,
               texname = '\\text{Mf}',
               lhablock = 'MASS',
               lhacode = [ 20 ])

Wh = Parameter(name = 'Wh',
               nature = 'external',
               type = 'real',
               value = 1,
               texname = '\\text{Wh}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

Wf = Parameter(name = 'Wf',
               nature = 'external',
               type = 'real',
               value = 1,
               texname = '\\text{Wf}',
               lhablock = 'DECAY',
               lhacode = [ 20 ])

WA = Parameter(name = 'WA',
               nature = 'external',
               type = 'real',
               value = 1,
               texname = '\\text{WA}',
               lhablock = 'DECAY',
               lhacode = [ 22 ])

EL = Parameter(name = 'EL',
               nature = 'internal',
               type = 'real',
               value = '1',
               texname = '\\text{EL}')

Yl = Parameter(name = 'Yl',
               nature = 'internal',
               type = 'real',
               value = '1',
               texname = '\\text{Yl}')

Yf = Parameter(name = 'Yf',
               nature = 'internal',
               type = 'real',
               value = '1',
               texname = '\\text{Yf}')

