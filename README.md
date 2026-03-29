# Earthquake Analysis and Catastrophe Modelling – Chennai

## Overview
This project analyses historical earthquake data and develops a simplified catastrophe modelling workflow for the Chennai region.

Due to the lack of granular exposure and ground motion data, a synthetic dataset was constructed to demonstrate the hazard–vulnerability–loss framework. The primary focus of the study is on exploratory data analysis (EDA) and understanding seismic patterns.

## Exploratory Data Analysis (EDA)

The USGS earthquake dataset (1976–2026) was analysed to understand key characteristics of seismic activity in India and surrounding regions.

### Key Findings

- Earthquake magnitudes follow the Gutenberg–Richter law, with a linear relationship observed in log-scale frequency plots  
- A completeness threshold of M ≥ 4.2 was identified, below which events are probably under-recorded  
- Most earthquakes are concentrated in the 4–5 magnitude range, with frequency declining for higher magnitudes  
- Significant clustering of events is observed along the Himalayan belt and Andaman–Sumatra region
- Earthquake depth and magnitude show low correlation, indicating independence  
- A spike in earthquake counts during 2004–2005 is attributed to aftershocks from the Sumatra–Andaman earthquake  

## Methodology

### 1. Simulation Framework
- Earthquake occurrences modelled using a Poisson process
- Magnitudes, depth, and locations generated using bootstrap sampling
- Synthetic dataset created to extend observations beyond historical limitations

### 2. Ground Motion Estimation
- Peak Ground Acceleration (PGA) estimated using an attenuation relationship
- Based on earthquake magnitude and hypocentral distance from Chennai

### 3. Exposure and Vulnerability
- Portfolio of buildings were simulated for Chennai
- Building types include Masonry and Reinforced Concrete (RC)
- Damage ratios assigned based on PGA levels

### 4. Loss Estimation
- Annual losses computed using simulated PGA and exposure
- Loss = Asset Value × Damage Ratio
- Loss exceedance curves constructed to assess risk

## Key Results

- Average Annual Loss (AAL): ₹9.7 million  
- 99th percentile loss: ₹69 million  
- 99.5th percentile loss: ₹142 million  

## Key Insights

- Most earthquakes produce negligible ground motion at the study location  
- Seismic risk is dominated by rare extreme events  
- Hazard and loss curves exhibit expected exponential decay behaviour  
- Data limitations can significantly impact tail risk estimation  

---

## Limitations

- Synthetic exposure data used due to lack of real portfolio data  
- Deterministic ground motion model (no variability)  
- Simplified vulnerability relationships  
- No variation in ground motion across assets  
- Limited number of extreme events

---

## Tools Used

- R (dplyr, ggplot2)
- USGS Earthquake Catalog

---

## Author

Sreyes K K
