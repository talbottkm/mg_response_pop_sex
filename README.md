1c_master_20250915.csv = all data for analyses; column headers below
band = individual identity
sex = male or female
primary_dose = first inoculation dose 750 CCU or 0 CCU of MG
secondary_dose = second inoculation dose 30, 100, 300, or 7000 CCU of MG
year = year of data collection 
population = population of origin (Arizona = AZ, Virginia = VA)
quantity## = number of MG copies by qPCR on ## days post initial inoculation
eyetot## = combined eye score on ## days post initial inoculation
prim_maxeye = maximum combined eyescore from the initial phase (7 to 41 days post initial inoculation)
prim_maxload = maximum MG load from the initial phase (7 to 41 days post initial inoculation)
prim_resist = load at 7 days post initial inoculation 
prim_resist10 = log10(prim_resist + 1)
prim_toleye = max combined eyescore between 7 & 14 days post initial inoculation, used to calculate tolerance
prim_inf = infection assignment during initial phase (1 = successfully infected, 0 = not successfully infected)
sec_maxeye = maximum combined eyescore from the second phase (46 to 63 days post initial inoculation)
sec_maxload = maximum MG load from the second phase (46 to 63 days post initial inoculation)
sec_resist = load at 49 days post initial inoculation 
sec_resist10 = log10(sec_resist + 1)
sec_toleye = max combined eyescore between 49 & 56 days post initial inoculation, used to calculate tolerance
sec_inf = infection assignment during second phase (1 = infected, 0 = not successfully infected)

1Cmanuscript_InitialAnalyses_20251120.R = R code for analyses investigating responses to initial inoculation, using birds receiving 0 or 750 CCU MG

1Cmanuscript_Main_20251120.R = R code for analyses for main manuscript, investigating responses to a second inoculation, using only birds that received 750 CCU MG during the initial inoculation
