# Input variables added by CED and their use

## Flexo phonon:
- `adcalc`: Set to 1 for adiabatic response calculation, 2 for diamagnetic (i.e., CRG) calculation (previously `userib`)
- `nogzero`: Remove G=0 term from electrostatic potential (previously `useria`)
- `joperloc`: Default (0) means joper uses local+nl; -1 mean joper uses only local  
- `pmpath`: 0 (default) is ICL path, 1 is PM
- `symfxe`: 1 to turn off symmetry in m_dfpt_looppert. For some reason needed for Flexo calculations

## Flexo metric:
- `metcalc`: Set to one to indicate that previous run was metric perturbation (previously `userid`)
-

## Velocity-force and naBEC:
- `vfstep`: Triggers velfrc, 1 for first pert, 2 for second (previously `useria`)
- `drudewt`: Turn on Drude weight, first digit is direction a, second is direction b (previously `useric`)
- `vlfrceta`: size of the imaginary part for causality (previously `userra`)

## Misc:


# Other notes:
- Command to push: `git push abinit-grp-develop flexo_cur_met_2:main`
- Command to fetch changes from gitlab: `git pull trunk develop`