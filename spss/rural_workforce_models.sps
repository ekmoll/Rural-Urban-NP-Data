* Encoding: UTF-8.
* =============================================================================.
* Rural Clinician Workforce Panel, 2018-2023: SPSS replication of the models.
*
* Reads the cleaned panel (data/AHRF_allphys_NP_PA_panel_2018_2023_clean.csv),
* reshapes it to one row per county-year-profession, and fits the Poisson models
* reported in the manuscript. Expected results from the Python scripts are given
* in comments under each model so the output can be compared line by line.
*
* HOW TO RUN
*   1. Change the path on the FILE HANDLE line below to where you saved the CSV.
*   2. Run All.
*
* NOTES ON MATCHING THE PYTHON RESULTS
*   - Python clusters standard errors by county. The SPSS equivalent is GEE
*     (REPEATED subject = county, independent working correlation, robust
*     covariance). Rate ratios should match exactly. Confidence limits can differ
*     slightly in the 3rd decimal because Python applies a small-sample
*     correction to the clustered errors and SPSS does not.
*   - Look at the Exp(B) column of "Parameter Estimates".
*   - Profession is coded 1 = NP, 2 = PA, 3 = Physician. SPSS uses the LAST
*     category as the reference, so physicians are the reference group.
*       rucc23                 = physician change per RUCC step (rate ratio)
*       [prof=1] * rucc23      = NP gradient relative to physicians
*       [prof=2] * rucc23      = PA gradient relative to physicians
*   - Year and region coefficients are adjustment terms only.
* =============================================================================.

FILE HANDLE panel /NAME='C:\path\to\AHRF_allphys_NP_PA_panel_2018_2023_clean.csv'.

GET DATA
  /TYPE=TXT
  /FILE=panel
  /ENCODING='UTF8'
  /DELIMITERS=","
  /QUALIFIER='"'
  /ARRANGEMENT=DELIMITED
  /FIRSTCASE=2
  /VARIABLES=
    fips A5
    county A40
    state A25
    region F1.0
    region_name A10
    rucc23 F1.0
    metro A5
    nonmetro A5
    rural A5
    rucc13 F1.0
    metro13 A5
    year F4.0
    population F10.0
    pop_imputed A5
    phys_all F8.0
    pcphys F8.0
    np F8.0
    pa F8.0
    apc F8.0
    total_three F8.0
    total_pc_capable F8.0
    phys_all_rate F12.6
    pcphys_rate F12.6
    np_rate F12.6
    pa_rate F12.6
    apc_rate F12.6
    phys_share F12.6
    np_share_all F12.6
    pa_share_all F12.6
    np_share F12.6
    apc_share F12.6.
EXECUTE.

* Expect 18,852 cases (3,142 counties x 6 years).
COMPUTE lnpop = LN(population).
COMPUTE is_nonmetro = (nonmetro = 'True').
COMPUTE year_c = year - 2018.
* Nonmetro urban population size (3 = 20,000+ is the reference) and adjacency.
RECODE rucc23 (4,5=3) (6,7=2) (8,9=1) (ELSE=SYSMIS) INTO usize.
VALUE LABELS usize 1 'Urban pop < 5,000' 2 'Urban pop 5,000-20,000' 3 'Urban pop 20,000+ (ref)'.
COMPUTE nonadj = ANY(rucc23, 5, 7, 9).
* RUCC as a category with RUCC 1 as the reference (recoded to 10 so it sorts last).
RECODE rucc23 (1=10) (ELSE=COPY) INTO rucc_cat.
VALUE LABELS rucc_cat 10 'RUCC 1 (ref)'.
EXECUTE.

DATASET NAME wide.
* Keep an untouched copy of the wide file for the primary care models at the end.
DATASET COPY wide_pcp.
DATASET ACTIVATE wide.

* -----------------------------------------------------------------------------.
* Long file, all physicians as comparator: 56,556 rows.
* -----------------------------------------------------------------------------.
VARSTOCASES
  /MAKE count FROM np pa phys_all
  /INDEX=prof(3)
  /KEEP=fips year year_c region rucc23 rucc_cat usize nonadj is_nonmetro lnpop.
VALUE LABELS prof 1 'NP' 2 'PA' 3 'Physician (ref)'.

* MODEL 1 (main model, Table 3).
* Expected (Python, county-clustered):
*   rucc23 Exp(B)           0.812 (0.793-0.832)   physicians -18.8% per step
*   [prof=1]*rucc23 Exp(B)  1.190 (1.165-1.215)   NPs  -3.4% per step, p < .001
*   [prof=2]*rucc23 Exp(B)  1.150 (1.126-1.176)   PAs  -6.6% per step, p < .001
*   NP rate ratio 0.966 (0.951-0.982); PA rate ratio 0.934 (0.917-0.951)
*   (= rucc23 Exp(B) x interaction Exp(B)).
GENLIN count BY prof year region WITH rucc23
  /MODEL prof rucc23 prof*rucc23 year region INTERCEPT=YES
    DISTRIBUTION=POISSON LINK=LOG OFFSET=lnpop
  /REPEATED SUBJECT=fips WITHINSUBJECT=year*prof SORT=YES
    CORRTYPE=INDEPENDENT ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO FIT SUMMARY SOLUTION (EXPONENTIATED).

* MODEL 2 (Supplemental Table S1): negative binomial, dispersion estimated.
* GEE cannot estimate the dispersion, so this runs without the REPEATED
* subcommand (model-based SEs). Compare rate ratios, not CIs.
* Expected rate ratios (Python NB2):
*   rucc23 0.852; [prof=1]*rucc23 1.146; [prof=2]*rucc23 1.129; alpha 0.518.
GENLIN count BY prof year region WITH rucc23
  /MODEL prof rucc23 prof*rucc23 year region INTERCEPT=YES
    DISTRIBUTION=NEGBIN(MLE) LINK=LOG OFFSET=lnpop
  /PRINT MODELINFO FIT SUMMARY SOLUTION (EXPONENTIATED).

* MODEL 3 (Supplemental Table S1): RUCC as categorical, RUCC 1 = reference.
* Expected physician rate ratios vs RUCC 1 (the rucc_cat terms):
*   RUCC 2 0.817, 3 0.735, 4 0.383, 5 0.669, 6 0.295, 7 0.420, 8 0.166, 9 0.187.
* NP and PA ratios = rucc_cat Exp(B) x matching [prof]*[rucc_cat] Exp(B):
*   NPs: 1.059, 1.121, 0.812, 1.200, 0.736, 0.995, 0.601, 0.729
*   PAs: 1.125, 1.070, 0.749, 1.004, 0.581, 0.773, 0.439, 0.603.
GENLIN count BY prof rucc_cat year region
  /MODEL prof rucc_cat prof*rucc_cat year region INTERCEPT=YES
    DISTRIBUTION=POISSON LINK=LOG OFFSET=lnpop
  /REPEATED SUBJECT=fips WITHINSUBJECT=year*prof SORT=YES
    CORRTYPE=INDEPENDENT ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO FIT SUMMARY SOLUTION (EXPONENTIATED).

* MODEL 4: did the profession gradients change over time?
* The year_c main effect is omitted because the year fixed effects absorb it.
* Expected: [prof=1]*rucc23*year_c Exp(B) 1.0002, p = .69
*           [prof=2]*rucc23*year_c Exp(B) 0.9974, p < .001.
GENLIN count BY prof year region WITH rucc23 year_c
  /MODEL prof rucc23 prof*rucc23 year region
    rucc23*year_c prof*year_c prof*rucc23*year_c INTERCEPT=YES
    DISTRIBUTION=POISSON LINK=LOG OFFSET=lnpop
  /REPEATED SUBJECT=fips WITHINSUBJECT=year*prof SORT=YES
    CORRTYPE=INDEPENDENT ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO FIT SUMMARY SOLUTION (EXPONENTIATED).

* MODEL 5 (Supplemental Table S2): nonmetro counties only (1,957 counties).
* Expected physician rate ratios: usize=1 (<5,000) 0.372 (0.326-0.424);
*   usize=2 (5,000-20,000) 0.729 (0.625-0.851); nonadj 1.534 (1.298-1.813).
* NPs (main x interaction): <5,000 0.708; 5,000-20,000 0.900; nonadj 1.377.
* PAs: <5,000 0.615; 5,000-20,000 0.789; nonadj 1.352.
* Interaction p-values: [prof=1]*nonadj p = .094; [prof=2]*nonadj p = .033.
TEMPORARY.
SELECT IF is_nonmetro = 1.
GENLIN count BY prof usize year region WITH nonadj
  /MODEL prof usize prof*usize nonadj prof*nonadj year region INTERCEPT=YES
    DISTRIBUTION=POISSON LINK=LOG OFFSET=lnpop
  /REPEATED SUBJECT=fips WITHINSUBJECT=year*prof SORT=YES
    CORRTYPE=INDEPENDENT ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO FIT SUMMARY SOLUTION (EXPONENTIATED).

* -----------------------------------------------------------------------------.
* Long file, PRIMARY CARE physicians as comparator (sensitivity analysis).
* -----------------------------------------------------------------------------.
DATASET ACTIVATE wide_pcp.
VARSTOCASES
  /MAKE count FROM np pa pcphys
  /INDEX=prof(3)
  /KEEP=fips year region rucc23 lnpop.
VALUE LABELS prof 1 'NP' 2 'PA' 3 'Primary care physician (ref)'.

* MODEL 6: primary care physician comparator.
* Expected: rucc23 Exp(B) 0.914 (0.903-0.925), PCPs -8.6% per step
*   [prof=1]*rucc23 1.055 (1.037-1.073), p < .001   (NPs -3.6% per step)
*   [prof=2]*rucc23 1.020 (1.000-1.040), p = .051   (PAs -6.8% per step).
GENLIN count BY prof year region WITH rucc23
  /MODEL prof rucc23 prof*rucc23 year region INTERCEPT=YES
    DISTRIBUTION=POISSON LINK=LOG OFFSET=lnpop
  /REPEATED SUBJECT=fips WITHINSUBJECT=year*prof SORT=YES
    CORRTYPE=INDEPENDENT ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO FIT SUMMARY SOLUTION (EXPONENTIATED).
