#1/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N SImpSamp

Enthought/Canopy_64bit/User/bin/python "ToricCodeProject2/bin/SingleShotLinearLengthSweep.py" $SGE_TASK_ID
