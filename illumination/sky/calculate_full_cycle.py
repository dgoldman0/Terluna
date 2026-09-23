from solver import *
from render_data import build_atlas
if __name__=='__main__':
 for a in [EARTH,MOON,MOON_ZERO]:
  f=solve(a,'standard',20.,85)
  build_atlas(f)
