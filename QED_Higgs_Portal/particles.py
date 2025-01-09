# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Mon 30 Dec 2024 02:28:48


from __future__ import division
from object_library import all_particles, Particle
import parameters as Param

import propagators as Prop

H = Particle(pdg_code = 25,
             name = 'H',
             antiname = 'H',
             spin = 1,
             color = 1,
             mass = Param.MH,
             width = Param.WH,
             texname = 'H',
             antitexname = 'H',
             charge = 0)

e = Particle(pdg_code = 11,
             name = 'e',
             antiname = 'e~',
             spin = 2,
             color = 1,
             mass = Param.Me,
             width = Param.ZERO,
             texname = 'e',
             antitexname = 'e~',
             charge = 1)

e__tilde__ = e.anti()

m = Particle(pdg_code = 13,
             name = 'm',
             antiname = 'm~',
             spin = 2,
             color = 1,
             mass = Param.MM,
             width = Param.ZERO,
             texname = 'm',
             antitexname = 'm~',
             charge = 1)

m__tilde__ = m.anti()

tt = Particle(pdg_code = 15,
              name = 'tt',
              antiname = 'tt~',
              spin = 2,
              color = 1,
              mass = Param.MTA,
              width = Param.ZERO,
              texname = 'tt',
              antitexname = 'tt~',
              charge = 1)

tt__tilde__ = tt.anti()

A = Particle(pdg_code = 22,
             name = 'A',
             antiname = 'A',
             spin = 3,
             color = 1,
             mass = Param.ZERO,
             width = Param.ZERO,
             texname = 'A',
             antitexname = 'A',
             charge = 0)

