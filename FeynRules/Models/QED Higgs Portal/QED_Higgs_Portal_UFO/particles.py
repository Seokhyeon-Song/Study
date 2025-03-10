# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Thu 26 Dec 2024 16:01:58


from __future__ import division
from object_library import all_particles, Particle
import parameters as Param

import propagators as Prop

h = Particle(pdg_code = 25,
             name = 'h',
             antiname = 'h',
             spin = 1,
             color = 1,
             mass = Param.Mh,
             width = Param.Wh,
             texname = 'h',
             antitexname = 'h',
             charge = 0)

chi = Particle(pdg_code = 38,
               name = 'chi',
               antiname = 'chi~',
               spin = 1,
               color = 1,
               mass = Param.Mchi,
               width = Param.ZERO,
               texname = 'chi',
               antitexname = 'chi~',
               charge = 0)

chi__tilde__ = chi.anti()

ele = Particle(pdg_code = 11,
               name = 'ele',
               antiname = 'ele~',
               spin = 2,
               color = 1,
               mass = Param.ME,
               width = Param.ZERO,
               texname = 'ele',
               antitexname = 'ele~',
               charge = 1)

ele__tilde__ = ele.anti()

mu = Particle(pdg_code = 13,
              name = 'mu',
              antiname = 'mu~',
              spin = 2,
              color = 1,
              mass = Param.MM,
              width = Param.ZERO,
              texname = 'mu',
              antitexname = 'mu~',
              charge = 1)

mu__tilde__ = mu.anti()

ta = Particle(pdg_code = 15,
              name = 'ta',
              antiname = 'ta~',
              spin = 2,
              color = 1,
              mass = Param.ML,
              width = Param.ZERO,
              texname = 'ta',
              antitexname = 'ta~',
              charge = 1)

ta__tilde__ = ta.anti()

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

