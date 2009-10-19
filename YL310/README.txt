This is a readme file. 

Scripts will be documented with examples and comments once exams are over.
You can always see the command line options with the -h switch.

=== For now, here are some very basic examples with terse comments: ===

-=> python gaussian_data_latest.py -r *.log 
(renames files systematically)

-=> python gaussian_data_latest.py -w *.g03 
(uploads to mysql - edit script with db name first)

-=> python mysql_thermo_data_latest.py -r -i 213{3,4,5,6,7,8,9} 
(runs thermo calcs and prints to screen)

-=> python mysql_thermo_data_latest.py -w -i 213{3,4,5,6,7,8,9} 
(runs thermo calcs and uploads to db)

-=> python enthalpy_calculation_latest.py -i 2133 
(calculates formation enthalpy, specify -w to upload to db)

-=> python spin_state_calc_latest.py -w *.mol 
(calculates spin states for all mol files in dir)

-=> python mysql_data_copy.py --from=YL310thermoname_B3LYP --field=symmetrynumber --to=YL310thermoname_B971 -i 21{4,5}{0,1,3} 
(copies symmetry numbers between databases for matching SMILES names and spin states only, for specified IDs in source db)

-=> python data_retrieval.py --csv=b3lyp_2 -i 213{4,5,6,7,8,9}
(retrieves data and dumps to CSV. Currently does the same table data as seen in my report, but includes number of Al, Ti, O, Cl molecules to facilitate sorting first)


=== Typical workflow: ===

1. python gaussian_data_latest.py -r *.log (renames files first so that mol filenames in db are nice)
2. python gaussian_data_latest.py -w *.g03 
3. python mysql_thermo_data_latest.py -w -i 213{3,4,5,6,7,8,9} 
4. python enthalpy_calculation_latest.py -i 213{3,4,5,6,7,8,9} 
