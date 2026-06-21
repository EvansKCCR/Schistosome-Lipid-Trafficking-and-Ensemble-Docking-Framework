# Post-MD receptor-ligand complex analysis report

Input workbook: `/mnt/data/CEH_trajectory_data.xlsx`
Mode: `ceh`
Summary workbook: `post_MD_complex_analysis_summary.xlsx`
Figures directory: `figures/`

## Outputs

- `post_MD_complex_analysis_summary.xlsx`: summary tables and long-form analysis outputs.
- `figures/`: comparative RMSD, ligand RMSD, Rg, SASA, RMSF, contact, DSSP, and CEH catalytic plots when available.

## Comparative index preview

| System     | Ligand     |   C_alpha_RMSD |   Lig_RMSD |    rGyr |   Mean_catalytic_score_0_to_5 |   Mean_catalytic_score_fraction |   Fully_productive_ok_percent |   Partial_catalytic_ok_percent |
|:-----------|:-----------|---------------:|-----------:|--------:|------------------------------:|--------------------------------:|------------------------------:|-------------------------------:|
| CEH_CE18_1 | CEH–CE18:1 |       0.158508 |   0.467407 | 2.51805 |                      0.346614 |                       0.0693227 |                             0 |                              0 |
| CEH_CE18_2 | CEH–CE18:2 |       0.185039 |   0.507587 | 2.53375 |                      0.063745 |                       0.012749  |                             0 |                              0 |
| CEH_CE20_4 | CEH–CE20:4 |       0.196993 |   0.418581 | 2.5348  |                      0.87251  |                       0.174502  |                             0 |                              0 |
| CEH_CE20_5 | CEH–CE20:5 |       0.150347 |   0.348223 | 2.52423 |                      1.72112  |                       0.344223  |                             0 |                              0 |

## CEH catalytic criteria preview

| System     | Ligand     |   Analysis_Start_ns |   N_frames |   Mean_catalytic_score_0_to_5 |   Mean_catalytic_score_fraction |   ser207_his508_ok_percent |   his508_glu388_oe2_ok_percent |   ser207_carbonylC_distance_ok_percent |   gly128N_O2_distance_ok_percent |   gly129N_O2_distance_ok_percent |   Triad_contact_ok_percent |   Attack_distance_ok_percent |   Oxyanion_any_ok_percent |   Oxyanion_both_ok_percent |   Partial_catalytic_ok_percent |   Fully_productive_ok_percent |
|:-----------|:-----------|--------------------:|-----------:|------------------------------:|--------------------------------:|---------------------------:|-------------------------------:|---------------------------------------:|---------------------------------:|---------------------------------:|---------------------------:|-----------------------------:|--------------------------:|---------------------------:|-------------------------------:|------------------------------:|
| CEH_CE18_1 | CEH–CE18:1 |                   0 |        251 |                      0.346614 |                       0.0693227 |                   2.39044  |                      32.2709   |                                0       |                                0 |                           0      |                   0.796813 |                      0       |                    0      |                          0 |                              0 |                             0 |
| CEH_CE18_2 | CEH–CE18:2 |                   0 |        251 |                      0.063745 |                       0.012749  |                   0.796813 |                       0        |                                5.57769 |                                0 |                           0      |                   0        |                      5.57769 |                    0      |                          0 |                              0 |                             0 |
| CEH_CE20_4 | CEH–CE20:4 |                   0 |        251 |                      0.87251  |                       0.174502  |                   0        |                       0        |                               87.251   |                                0 |                           0      |                   0        |                     87.251   |                    0      |                          0 |                              0 |                             0 |
| CEH_CE20_5 | CEH–CE20:5 |                   0 |        251 |                      1.72112  |                       0.344223  |                 100        |                       0.398406 |                               36.6534  |                                0 |                          35.0598 |                   0.398406 |                     36.6534  |                   35.0598 |                          0 |                              0 |                             0 |