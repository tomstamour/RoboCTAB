# RoboCTAB - Automated DNA Extraction Protocol

An automated CTAB (Cetyltrimethylammonium Bromide) DNA extraction protocol designed for Opentrons liquid handling robots. This protocol automates the traditional manual CTAB DNA extraction method, adapted from Doyle & Doyle (1987), enabling high-throughput and reproducible DNA extraction from plant tissues.

**Read the full method article here:** https://www.mdpi.com/2223-7747/14/15/2263

**Please cite this automated method as follow:**
***"Recent advances in robotic DNA extraction have improved throughput for large-scale genotyping studies (Boucher St-Amour et al., 2025)."***

## Key Features

- **Multi-plate Processing**: Handle 1 to 4 sample plates (96-well format) in a single run
- **Flexible Mixing Options**: Choose between pipette mixing, bubble mixing, or no mixing for chloroform step
- **Automatic Reagent Calculation**: Built-in calculations for all reagents based on sample number
- **Configurable Parameters**: Easy customization of volumes, labware, and processing options
- **Time Estimation**: Automatic calculation of protocol duration for key steps

## Requirements

- Opentrons Robot
- P300 Multi-Channel Pipette
- 1.2ml Deep-well plates (96-well format) for samples
- 300µL tip racks (filtered or non-filtered)
- Reservoirs (minimum 200ml capacity)
- Centrifuge capable of 6000 rpm
- INcubator or water bath set to 65°C

## Configuration

1. Download the `RoboCTAB.py` file to your Opentrons protocols folder
2. Before running the protocol, modify the configuration parameters in lines 6-23 of `RoboCTAB.py`
3. Load the protocol in the Opentrons App
4. Verify all labware definitions are recognized

### Protocol Execution

The protocol consists of several main steps:

1. **TE Buffer Addition** - Initial buffer dispensing to samples
2. **Tissue Grinding** - Manual grinding using tissue lyser
3. **Extraction Buffer Addition** - Hot extraction buffer dispensing
4. **Incubation** - 30 to 60 minutes at 65°C
5. **Chloroform Extraction** - Organic phase separation with mixing
6. **Centrifugation** - Phase separation 
7. **Supernatant Transfer** - Aqueous phase collection
8. **Isopropanol Precipitation** - DNA precipitation
9. **Isopropanol Wash** - Pellet washing
10. **Ethanol Wash** - Final washing step
11. **Drying** - Ethanol evaporation
12. **Elution** - DNA resuspension in elution buffer

## Citation

This automated protocol is based on the CTAB DNA extraction method originally described by:

Doyle, J.J. and Doyle, J.L. (1987). A rapid DNA isolation procedure for small quantities of fresh leaf tissue. *Phytochemical Bulletin*, 19, 11-15.

## Contact

For questions, issues, or contributions, please contact:
tomstamour@hotmail.com

---

*This protocol was developed for research purposes. Users are responsible for validation and optimization for their specific applications.
