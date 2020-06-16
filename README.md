# MATLAB Code for "Actigraphy-based Parameter Tuning Process for Adaptive Notch Filter and Circadian Phase Shift Estimation" in Chronobiology International
Matlab codes for the adaptive notch filter (ANF) algorithm and ANF parameter tuning algorithm proposed in Chronobiology International paper "Actigraphy-based Parameter Tuning Process for Adaptive Notch Filter and Circadian Phase Shift Estimation".
# MATLAB Download and Installation
The latest release MATLAB R2020a can be downloaded here:
https://www.mathworks.com/products/get-matlab.html?s_tid=gn_getml

The installation of MATLAB could follow the instruction in MathWorks:
https://www.mathworks.com/help/install/install-products.html

# Invoking the MATLAB Codes (.m) and Simulink files (.mdl)
The five Simulink files "ANF_1st.mdl", "ANF_2nd.mdl", "ANF_3rd.mdl", "ANF_4th.mdl" and "ANF_5th.mdl" are used to run the 1st-order, 2nd-order to 5th-order ANF, respectively.

The MATLAB codes (.m) can invoke the Simulink files directly. The file "Adaptive_notch_filter.m" runs 1st-order, 2nd-order, 3rd-order, 4th-order, and 5th-order ANF one by one, plots and compares their estimated value and circadian phase.

The file "Evolutionary_Strategy_ANF_mutation.m" uses the Evolution Strategy to optimize the ANF parameters by minimizing the objective function proposed in our Chronobiology International paper. The final results of 'zeta', 'gamma_omg', and 'gamma_d' give the optimal damping factor, adaptation rates of notch frequency and constant bias, respectively. The algorithm of Evolution Strategy can be found in:
https://www.sft.asso.fr/Local/sft/dir/user-3775/documents/actes/journeessft/metti_5_2012/Lectures&Tutorials-Texts/Text-T2-Ruffio.pdf
