from distutils.core import setup

setup(name='whisker_sim',
      version='0.1',
      description='simulates whiskers using trep',
      author='Elliot Hevel',
      author_email='elliothevel2013@u.northwestern.edu',
      package_dir={'whisker_sim': 'src'},
      packages=['whisker_sim', 'whisker_sim.collisions',
          'whisker_sim.simulation',
          'whisker_sim.whisker','whisker_sim.collisions.surfaces'],
      package_data={'':['simulation/rathead.stl']})

