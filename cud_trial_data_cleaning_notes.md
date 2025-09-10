# CUD Trial Data Cleaning Notes

## Timepoints (`Time` Variable)

| Timepoint Name             | Code |
|----------------------------|------|
| Prescreening Assessment 1  | 1    |
| Prescreening Assessment 2  | 2    |
| Prescreening Assessment 3  | 3    |
| Prescreening Assessment 4  | 4    |
| Medical Screening          | 5    |
| Preparation 1              | 6    |
| Preparation 2              | 7    |
| Preparation 3              | 8    |
| Preparation 4              | 9    |
| Preparation 5              | 10   |
| MRI 1                      | 11   |
| Drug Administration        | 12   |
| Integration / MRI2         | 13   |
| Follow-up 1                | 14   |
| Follow-up 2                | 15   |
| Follow-up 3                | 16   |
| Follow-up 4                | 17   |
| 90-day Assessment          | 18   |
| 180-day Assessment         | 19   |
| Blood Pressure Check       | 20   |

---

## Main Analysis Outcome (MMRM)

| `Time_Recoded` | Time Values       | Period        | Peter-Created Variable     |
|----------------|-------------------|---------------|-----------------------------|
| 1              | 1                 | Baseline      | `AB_Per_PreContact`        |
| 2              | 2–6               | Screening     | `AB_Per_PreTX`             |
| 3              | 7–12              | Preparation   | `AB_Per_Preps`             |
| 4              | 13–17             | Integration   | `AB_Per_PostDrugTX`        |
| 5              | 18                | 90 Day        | `AB_Per_Follow_Ups`        |
| 6              | 19                | 180 Day       | `AB_Per_Follow_Ups`        |

---

## Secondary Analysis Outcome – Kaplan-Meier + Cox

| Days Variable | Recorded Timepoint | Represented Period                         |
|---------------|--------------------|--------------------------------------------|
| Prescreen     | Prescreen 1        | Baseline                                   |
| Preparation   | Prep Session 1     | Screening                                  |
| Drug          | Administration     | Preparation                                |
| Integration   | Integration 1?     | Full Post-Drug Period (Integration + F/U)? |

---

## Observation Variables (Description)

>
> 1. Number of days under observation during the pre-screen period  
> 2. Number of days under observation during the preparation period  
> 3. Number of days between the last preparation period and the drug administration  
> 4. Number of days under observation during the integration period”

---

## Issues with Survival Data

There are 6 participants missing `DAYS` values:

- **Placebo**: 1012, 1037  
- **Psilocybin**: 1011, 1021, 1024, 1025, 1039

*1037* is accounted for. 
The remaining 6 were not accounted for and were handled using the `Integration` variable instead in the survival analysis dataset.

---

