# **Hybrid Precoding of mmWave Communication**

---

## **Abstract**

Millimeter Wave (mm-Wave) communication can significantly enhance 5G cellular technology's capacity to meet increasing demands of high data rate software application and media access. However, its small wavelength reduces the propagation range and is sensitive to obstacles, causing notable pathloss.

To address the high path loss of mm-Wave communication, traditional precoding methods such as fully digital precoding are cost-inefficient due to high power consumption from dedicated analog-to-digital converters (ADCs) per antenna. Hybrid precoding is introduced as a solution to reduce RF chains by mapping multiple antennas to a single RF chain. It combines the merits of analog and digital precoding: analog phase shifters steer beams directionally, while digital baseband precoding allows flexible phase assignment and power allocation.

This project assesses hybrid precoding algorithms for optimal performance. The hybrid precoding structure is reformulated as a matrix decomposition optimization problem using three methods: Optimal Unconstrained Precoding, Beam Steering, and Orthogonal Matching Pursuit (OMP). Their performance is compared in terms of spectral efficiency, demonstrating the effectiveness of the OMP-based hybrid precoding approach.

---

## **1. Background of mmWave**

### 1.1 What is mmWave?
- Frequency range: 30 GHz – 300 GHz
- Wavelength: 1 mm – 10 mm
- Enables large-scale antenna arrays and wide bandwidth (~250 GHz)
- Supports high data-rate communication (Gbps)
- In Hong Kong, current 5G bands include sub-6GHz and mmWave (e.g., n257, n261)

### 1.2 Problems of mmWave
- **Poor penetration**: high signal attenuation
- **Sensitive to obstacles**: scattering effects (reflection, refraction, diffraction)
- **Multipath propagation**: leads to spatial sparsity in channel paths

### 1.3 Beamforming
- Uses multiple antennas to focus signal energy in specific spatial directions
- Higher directionality with more antennas
- Enhances signal strength and mitigates interference

### 1.4 Existing Solutions
1. **Analog Beamforming**  
   - Uses phase shifters and RF chains
   - Low hardware complexity, but supports only single-stream transmission
2. **Fully Digital Precoding**  
   - Uses dedicated ADC per antenna
   - High performance but high power consumption and cost

### 1.5 Proposed Solution: Hybrid Precoding
- Combines analog and digital processing
- Reduces RF chains while maintaining performance
- Balances hardware complexity and spectral efficiency

---

## **2. System Model of Hybrid Precoding**

### 2.1 System Overview
- Single transmitter and single receiver
- Hybrid precoding splits processing into:
  - **Digital baseband precoder** (\( \mathbf{F}_{BB} \))
  - **Analog RF precoder** (\( \mathbf{F}_{RF} \))

### 2.2 Transmitter Architecture
\[
\mathcal{X} = \mathbf{F}_{RF} \mathbf{F}_{BB} \mathbf{S}
\]
where:
- \( \mathbf{S} \): data symbol vector
- \( \mathbf{F}_{BB} \): baseband precoder
- \( \mathbf{F}_{RF} \): analog precoder

### 2.3 Receiver Architecture
\[
\tilde{\mathbf{y}} = \sqrt{\rho} \mathbf{W}_{BB}^* \mathbf{W}_{RF}^* \mathbf{H} \mathbf{F}_{RF} \mathbf{F}_{BB} \mathbf{S} + \mathbf{W}_{BB}^* \mathbf{W}_{RF}^* \mathbf{n}
\]

### 2.4 Channel Model
Extended Saleh-Valenzuela cluster model:
\[
\mathbf{H} = \sqrt{\frac{N_t N_r}{N_{cl} N_{ray}}} \sum_{i=1}^{N_{cl}} \sum_{\ell=1}^{N_{ray}} \alpha_{i\ell} \mathbf{a}_r(\theta_{i\ell}^r) \mathbf{a}_t(\theta_{i\ell}^t)^H
\]

### 2.5 Antenna Array
- Uniform Linear Array (ULA) used
- Array response vector:
\[
\mathbf{a}(\theta) = \frac{1}{\sqrt{N}} \left[ 1, e^{j k d \sin \theta}, \dots, e^{j (N-1) k d \sin \theta} \right]^T
\]

---

## **3. Algorithm Design**

### 3.1 Optimal Unconstrained Precoding
- Uses Singular Value Decomposition (SVD) of channel matrix \( \mathbf{H} \)
- Selects dominant eigenvectors for precoding:
\[
\mathbf{F}_{opt} = \mathbf{V}_{1:N_s}, \quad \mathbf{W}_{opt} = \mathbf{U}_{1:N_s}
\]

### 3.2 Hybrid Precoding with Orthogonal Matching Pursuit (OMP)
**Objective:**  
\[
\min \| \mathbf{F}_{opt} - \mathbf{F}_{RF} \mathbf{F}_{BB} \|_F
\]
subject to analog constraints.

**OMP Steps:**
1. Initialize \( \mathbf{F}_{RF} \) as empty, residual \( \mathbf{F}_{res} = \mathbf{F}_{opt} \)
2. Select array response vector maximizing correlation with residual
3. Update \( \mathbf{F}_{RF} \) and compute \( \mathbf{F}_{BB} \) via least squares
4. Update residual and iterate

### 3.3 Beam Steering
- Brute-force dictionary search over beam directions
- Maximizes beam gain:
\[
(w^{opt}, f^{opt}) = \arg \max |w^* H f|^2
\]

---

## **4. Performance Analysis**

### 4.1 Simulation Parameters
| Parameter          | Value                    |
|---------------------|--------------------------|
| Transmit antennas   | 256                      |
| Receive antennas    | 64                       |
| RF chains           | 4–10                     |
| Data streams        | 2                        |
| SNR range           | -40 dB to 40 dB          |
| Channel clusters    | 8                        |
| Paths per cluster   | 10                       |

### 4.2 Results

**Global Spectral Efficiency Trends:**
- OMP-based hybrid precoding approaches optimal SVD performance
- Beam steering shows lower efficiency, especially at high SNR

**Hardware Efficiency:**
- With ≥8 RF chains, OMP achieves near-optimal performance
- Significant reduction in RF chains compared to fully digital (256 → 8)

**Beamforming Efficiency:**
- 20% spectral efficiency improvement with narrow, finely steered beams
- More antennas lead to better spatial resolution and interference suppression

---

## **5. Conclusion**

Hybrid precoding with OMP offers:
1. **Near-optimal performance** with significantly lower hardware cost
2. **Reduced RF chains** (e.g., 8 vs. 256 in fully digital)
3. **Efficient beamforming** through combined analog and digital processing
4. **Improved spectral efficiency** compared to analog beam steering

This work demonstrates that hybrid precoding is a practical and efficient solution for mmWave MIMO systems, balancing performance, cost, and complexity.

---

## **References**

1. O. E. Ayach et al., “Spatially sparse precoding in millimeter wave MIMO systems,” *IEEE Trans. Wireless Commun.*, 2014.
2. O. E. Ayach et al., “The capacity optimality of beam steering in large millimeter wave MIMO systems,” *IEEE WCNC*, 2012.
3. X. Yu et al., “Hybrid precoding design in millimeter wave MIMO systems: An alternating minimization approach,” *IEEE GLOBECOM*, 2015.
4. CMHK, “5G Network Trial Report,” OFCA, 2020.

---

**Thank you!**

--- 
