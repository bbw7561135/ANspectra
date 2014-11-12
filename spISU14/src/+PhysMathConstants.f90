!**********************************************************************!
MODULE PhysMathConstants
!----------------------------------------------------------------------!
! This module contains several mathematical and a lot of physical
! constants
!
! REFERENCES
!
! [1] T.G. Trippe et al., Rev. Mod. Phys. 48 (1976) S1-S245
! [2] R.L. Kelly et al.,  Rev. Mod. Phys. 52 (1980) S1-S286
! [3] C.G. Wohl et al.,  Rev. Mod. Phys. 56 (1984) S1-S299
! [4] D.E. Groom et al., Eur. Phys. Jour. C 15 (2000) 1-878
! [5] K. Hagiwara et al., Phys. Rev. D 66 (2002) 010001
! [6] S. Eidelman et al., Phys. Lett. B 592 (2004) 1-1109
! [7] W.-M. Yao et al., J. Phys. G. (Nucl. Part. Phys.) 33 (2006) 1-1232
! [8] J. Phys. G. (Nucl. Part. Phys.) 37 (2010) 075021
!----------------------------------------------------------------------!
!                                               ITEP, Moscow 2008/05/24!
!                                               JINR, Dubna  2014/10/15!
!**********************************************************************!
real,parameter:: &
!----------------------------------------------------------------------!
!NUMERICAL CONSTANTS
!----------------------------------------------------------------------!
    pi       = 3.141592653589793238462643383279503,&
    deg      = 180/pi,&
    dtr      = pi/180,&                                                !degrees to radians coefficient

    Infty    = 1d+10,&                                                 !"phisical infty"
    Precision= 0,&!epsilon(1.0),
!----------------------------------------------------------------------!
!PARTICLES MASSES [GeV/c^2]                                            !ATTENTION! Values are very old or just odd, some of them are very long away from PDG-2012!
!----------------------------------------------------------------------!
    m_e    = 5.1099892d-04,    mm_e    = m_e    **2,&                  !electron
    m_mu   = 1.05658369d-01,   mm_mu   = m_mu   **2,&                  !muon
    m_tau  = 1.77699,          mm_tau  = m_tau  **2,&                  !tau-lepton
    m_pi0  = 1.349766d-01,     mm_pi0  = m_pi0  **2,&                  !pi^0  light unflavored meson
    m_pi   = 1.3957018d-01,    mm_pi   = m_pi   **2,&                  !pi^+- light unflavored meson
    m_K    = 4.93677d-01,      mm_K    = m_K    **2,&                  !K^+-  strange          meson
    m_K0   = 4.97648d-01,      mm_K0   = m_K0   **2,&                  !K^0   strange          meson
    m_a1   = 1.230,            mm_a1   = m_a1   **2,&                  !a_1(1260)              meson
    m_Dpm  = 1.8694,           mm_Dpm  = m_Dpm  **2,&                  !D^+-  charmed          meson
    m_D0   = 1.8646,           mm_D0   = m_D0   **2,&                  !D^0   charmed          meson
    m_Dsm  = 1.9686,           mm_Dsm  = m_Dsm  **2,&                  !D^-_s charmed strange  meson
    m_Bpm  = 5.279,            mm_Bpm  = m_Bpm  **2,&                  !B^+-  bottom           meson
    m_B0   = 5.2794,           mm_B0   = m_B0   **2,&                  !B^0   bottom           meson
    m_rho  = 7.758d-01,        mm_rho  = m_rho  **2,&                  !rho                    meson
    m_r1   = 1.45,             mm_r1   = m_r1   **2,&                  !\rho(1450)             meson
    m_w0   = 7.84d-01,         mm_w0   = m_w0   **2,&                  !\omega(782)            meson
    m_w1   = 1.419,            mm_w1   = m_w1   **2,&                  !\omega(1420)           meson
    m_f0   = 1.019,            mm_f0   = m_f0   **2,&                  !\fhi(1020)             meson
    m_p    = 9.3827203d-01,    mm_p    = m_p    **2,&                  !proton
    m_n    = 9.3956536d-01,    mm_n    = m_n    **2,&                  !neutron
    m_I    = 0.5*(m_p+m_n),    mm_I    = m_I    **2,&                  !isoscalar nucleon
    m_L    = 1.115683,         mm_L    = m_L    **2,&                  !\Lambda             baryon
    m_Sp   = 1.18937,          mm_Sp   = m_Sp   **2,&                  !\Sigma^+            baryon
    m_Spp_c= 2.4525,           mm_Spp_c= m_Spp_c**2,&                  !\Sigma^++_c         baryon
    m_S0   = 1.192642,         mm_S0   = m_S0   **2,&                  !\Sigma^0            baryon
    m_Sm   = 1.197449,         mm_Sm   = m_Sm   **2,&                  !\Sigma^-            baryon
    m_Lp_c = 2.2849,           mm_Lp_c = m_Lp_c **2,&                  !\Lambda^+_c charmed baryon

    m_u    =( 1.5+  4.0)/2000, mm_u    = m_u    **2,&                  !u quark (3.250d-03)
    m_d    =( 4.0+  8.0)/2000, mm_d    = m_d    **2,&                  !d quark (6.000d-03)
    m_s    =(80.0+130.0)/2000, mm_s    = m_s    **2,&                  !s quark (1.050d-01)
    m_c    = 1.3,              mm_c    = m_c    **2,&                  !c quark
    m_b    =( 4.1+  4.4)*0.5,  mm_b    = m_b    **2,&                  !b quark (4.250d+00)
    m_t_DO = 1.743d+02,&                                               !t quark (Direct observation)
    m_t_SM = 1.781d+02,&                                               !t quark (SM electroweak fit)
    m_t    = m_t_DO,           mm_t    = m_t    **2,&                  !t quark (our current choice)

    m_W    = 8.0425d+01,       mm_W    = m_W    **2,&                  !W gauge boson
    m_Z    = 9.11876d+01,      mm_Z    = m_Z    **2,&                  !Z gauge boson
    m_D1232= 1.232,            mm_D1232= m_D1232**2,&                  !D(1232) resonance
!----------------------------------------------------------------------!
!CABIBBO-KOBAYASHI-MASKAWA QUARK MIXING MATRIX
!----------------------------------------------------------------------!
    Vl_ud= 9.739d-01, Vr_ud= 9.751d-01, Va_ud= 0.5*(Vl_ud+Vr_ud),&
    Vl_us= 2.210d-01, Vr_us= 2.270d-01, Va_us= 0.5*(Vl_us+Vr_us),&
    Vl_ub= 2.900d-03, Vr_ub= 4.500d-03, Va_ub= 0.5*(Vl_ub+Vr_ub),&
    Vl_cd= 2.110d-01, Vr_cd= 2.270d-01, Va_cd= 0.5*(Vl_cd+Vr_cd),&
    Vl_cs= 9.730d-01, Vr_cs= 9.744d-01, Va_cs= 0.5*(Vl_cs+Vr_cs),&
    Vl_sb= 3.900d-02, Vr_sb= 4.400d-02, Va_sb= 0.5*(Vl_sb+Vr_sb),&
    Vl_td= 4.800d-03, Vr_td= 1.400d-02, Va_td= 0.5*(Vl_td+Vr_td),&
    Vl_ts= 3.700d-02, Vr_ts= 4.300d-02, Va_ts= 0.5*(Vl_ts+Vr_ts),&
    Vl_tb= 9.990d-01, Vr_tb= 9.992d-01, Va_tb= 0.5*(Vl_tb+Vr_tb),&
!----------------------------------------------------------------------!
    G_D1232= 1.200d-01,&                                               !Full width of D(1232) resonance
    B_tm3  = 1.737d-01,&                                               !\tau_{\mu3} decay branching ratio
!----------------------------------------------------------------------!
!OTHER CONSTANTS
!----------------------------------------------------------------------!
    barn      = 1.000d-24,&                                            ![cm^2]
    millibarn = 1.000d-03*barn,&
    microbarn = 1.000d-06*barn,&
    nanobarn  = 1.000d-09*barn,&
    picobarn  = 1.000d-12*barn,&
    Fermi     = 1.000d-13,&                                            ![cm]

    alpha_e   = 1.0/1.37035999679d+02,&                                !Fine-structure constant
    N_Avogadro= 6.02214179d+23,&                                       !Avogadro constant [1/mol]
    C_Eiler   = 5.77215665d-01,&                                       !Eiler constant
    G_Fermi   = 1.16637d-05,&                                          !Fermi constant [(hc)^3/GeV^2]
    E_GRE     = 0.5*m_W**2/m_e,&                                       !Glashow resonance energy
    UAMU      = 931494.04380,&                                         !Unified Atomic Mass Unit [keV]
    AMU       = 1.66053873d-24,&                                       !Atomic unit mass [g]
    hbar      = 6.582118890d-22,&                                      !Planck constant, reduced [MeV s]
    hbarc     = 1.973269631d-01*Fermi,&                                ![GeV cm]
    ce        = 2.99792458d10,&                                        ![cm/s] !Stupid name. Never used
    hc2       = 3.893792920d-28,&                                      !Convertion constant [cm^2 GeV^2]
    lambda_e  = 3.8616d-11,&                                           !Compton vawelength of electron [cm]
    r_e       = 2.817940285d-13,&                                      !Classical electron radius [cm]
    lnR       = 1.890d+02,&                                            !Value of the radiation logarithm
    mu_p      = 2.792847337,&                                          !Proton  total magnetic moment
    mu_n      =-1.91304272,&                                           !Neutron total magnetic moment
    k_p       = mu_p-1.0,&                                             !Proton anomalous magnetic moment
    k_n       = mu_n,&                                                 !Neutron anomalous magnetic moment
    sin2W     = 2.3120d-01,&                                           !sin^2(Weinberg angle)
    cosC      = 9.7377d-01,&                                           !cos(Cabibbo angle)
    c2C       = cosC**2,&                                              !cos^2(Cabibbo angle)
    s2C       = 1.0-c2C,&                                              !sin^2(Cabibbo angle)
    R_Earth   = 6.371d+08,&                                            !Radius of the Earth in the PREM [cm]
    R_Sun     = 6.9598d+10,&                                           !Equatorial radius of the Sun [cm]
    E_GZK     = 5d+10,&                                                !GZK cuttof energy
    n_w       = 1.34!32986                                             !Refractive index of water
!----------------------------------------------------------------------!
real &
    NucParam(8,-1:18)/&
    1,   2, 0,      0,      0,      0,      0,    0,&                  !D2
    1,   1, 0,      0,      0,      0,      0,    0,&                  !H
    3,   7, 0.0151, 0.0151, 0.1690, 0.1690, 0.01, 0.01,&               !Li
    6,  12, 0.0256, 0.0257, 0.2210, 0.2210, 0.01, 0.01,&               !C
    8,  16, 0.0266, 0.0266, 0.2250, 0.2250, 0.01, 0.01,&               !O
    12, 24, 0.0289, 0.0290, 0.2350, 0.2350, 0.01, 0.01,&               !Mg
    14, 28, 0.0299, 0.0300, 0.2390, 0.2390, 0.01, 0.01,&               !Si
    20, 40, 0.0329, 0.0330, 0.2510, 0.2510, 0.01, 0.01,&               !Ca
    26, 56, 0.0361, 0.0330, 0.2630, 0.2510, 0.01, 0.01,&               !Fe
    28, 59, 0.0361, 0.0346, 0.2630, 0.2570, 0.01, 0.01,&               !Ni
    18, 40, 0.0350, 0.0307, 0.2590, 0.2420, 0.01, 0.01,&               !Ar
    50,119, 0.0391, 0.0300, 0.2740, 0.2390, 0.01, 0.01,&               !Sn
    82,207, 0.0417, 0.0315, 0.2830, 0.2450, 0.01, 0.01,&               !Pb
     7, 14, 0.0261, 0.0262, 0.2230, 0.2230, 0.01, 0.01,&               !N
     9, 19, 0.0283, 0.0284, 0.2325, 0.2325, 0.01, 0.01,&               !F
    10, 20, 0.0278, 0.0278, 0.2300, 0.2300, 0.01, 0.01,&               !Ne
    13, 27, 0.0294, 0.0295, 0.2390, 0.2390, 0.01, 0.01,&               !Al
    17, 36, 0.0314, 0.0315, 0.2450, 0.2450, 0.01, 0.01,&               !Cl
    29, 64, 0.0361, 0.0346, 0.2630, 0.2570, 0.01, 0.01,&               !Cu
    35, 80, 0.0381, 0.0315, 0.2703, 0.2450, 0.01, 0.01/                !Br
!    Z,  A,  E_bin[nu/an],  P_Fermi[nu/an],T_Fermi[nu/an]              !T_Fermi is random!!! P_Fermi is equal P_FeMAX!
!----------------------------------------------------------------------!
endMODULE PhysMathConstants
!**********************************************************************!
!                        ========================
!                        Z-A COMPOSITION OF WATER
!                        ========================
!     ------------- ------- -------- -------- --------- --------- 
!    | Experiment  |  <Z>  |  <A>   | <Z/A>  | <Z^2/A> | Density |
!    |  Location   |       |        |        |         | (g/cm^3)|
!     ------------- ------- -------- -------- --------- --------- 
!    | Suruga-bay, | 7.433 | 14.787 | 0.5530 |  3.760  |  1.027  |
!    | Shizuokaken |       |        |        |         |         |
!    |     [1]     |       |        |        |         |         |
!     ------------- ------- -------- -------- --------- --------- 
!    |  Black sea  | 7.430 | 14.780 | 0.5530 |  3.760  |    -    |
!    |     [2]     |       |        |        |         |         |
!     ------------- ------- -------- -------- --------- --------- 
!    | Carib. sea  |   -   |   -    |   -    |    -    |    -    |
!    |     [3]     |       |        |        |         |         |
!     ------------- ------- -------- -------- --------- --------- 
!    |    Hawaii   | 7.468 | 14.869 | 0.5525 |  3.770  |  1.027  |
!    |     [4]     |       |        |        |         |         |
!     ------------- ------- -------- -------- --------- --------- 
!    |    Pylos    |   -   |    -   |   -    |    -    |    -    |
!    |     [5]     |       |        |        |         |         |
!     ------------- ------- -------- -------- --------- --------- 
!    | lake Geneva |            pure water             |    1    |
!    |     [6]     |                                   |         |
!    |------------- ----------------------------------- ---------|
!    | lake Baikal |            pure water             |    1    |
!    |     [7]     |                                   |         |
!     ------------- ----------------------------------- --------- 
!                        ========================
!                        Z-A COMPOSITION OF ROCKS
!                        ========================
!     ------------- ------- -------- -------- --------- --------- 
!    | Laboratory  |  <Z>  |  <A>   | <Z/A>  | <Z^2/A> | Density |
!    | (substance) |       |        |        |         |(g/cm**3)|
!     ------------- ------- -------- -------- --------- --------- 
!    | Mont Blanc  | 10.34 | 20.94  | 0.4940 |  5.120  |  2.60   |
!    |     [8]     |       |        |(0.4156-| (3.667- |         |
!    |             |       |        | 0.4999)| 10.302) |         |
!     ------------- ------- -------- -------- --------- --------- 
!    |    Utah     |   -   |   -    |   -    |  5.650  |  2.47   |
!    |     [9]     |       |        |        | (5.540- | (2.40-  |
!    |             |       |        |        |  5.720) |  2.66)  |
!     ------------- ------- -------- -------- --------- --------- 
!    |    ERPM     |   -   |   -    | 0.4990 |  5.530  |  2.713  |
!    |    [10]     |       |        |+/-0.001| +/-0.01 |+/-0.007 |
!     ------------- ------- -------- -------- --------- --------- 
!    |    KGF      |  ~13  |  ~27   | 0.4950 |  6.310  |  3.050  |
!    |    [11]     |(12.8- | (26.9) |        |         | (2.98-  |
!    |             | 12.93)|        |        |         |  3.10)  |
!     ------------- ------- -------- -------- --------- --------- 
!    |   Frejus    |   -   |   -    | 0.5000 |  5.870  |  2.730  |
!    |    [12]     |       |        |        | +/-0.02 | +/-0.01 |
!     ------------- ------- -------- -------- --------- --------- 
!    | Gran Sasso  | 11.40 |   -    | 0.4980 |  5.700  |  2.710  |
!    |    [13]     |       |        |        |         | +/-0.05 |
!     ------------- ------- -------- -------- --------- --------- 
!    |   Baksan    | 11.90 |   -    | 0.4950 |  5.880  |  2.700  |
!    |    [14]     |       |        |        |         | +/-0.03 |
!     ------------- ------- -------- -------- --------- --------- 
!    |Standard rock|   11  |   22   | 0.5000 |  5.500  |  2.65   |
!    |    [15]     |       |        |        |         |         |
!     ------------- ------- -------- -------- --------- --------- 
!
! REFERENCES AND FOOTNOTES
!
!  [1] S. Higashi, T. Kitamura, S. Miyamoto et al., "Cosmic-Ray Inten-
!      sities Under Sea-Water at Depths Down to 1400 m.," Il Nuovo Ci-
!      mento, 43 A (1966) 334.
!  [2] Yu.N. Vavilov et al., Bull. Acad. of Sci. of the USSR 34 (1970)
!      1759.
!  [3] V.M. Fyodorov,  V.P. Pustovetov,  Yu.A. Trubkin and  A.V. Kiri-
!      lenkov, in Proc. of  the 19'th  International Cosmic Ray Confe-
!      rence,  La  Jolla,  California,  1985,  edited  by   F.C.Jones,
!      J.Adams and G.M. Mason  (NASA Conference  Publication No. 2376)
!      (Goddard  Space Flight Center, Greenbelt, MD,  Scientific  and 
!      Technical Information Branch, NASA, U.S. GPO, Washington, D.C.,
!      1985), Vol. 8, p. 39.
!  [4] J. Babson,  B. Barish,  R. Becker-Szendy  et  al.  (The  DUMAND
!      Collaboration),  Cosmic Ray Muons in the  Deep Ocean, HDC-1-89;
!      Available from Hawaii DUMAND Center,  University of Hawaii, Ho-
!      nolulu, 1989; Phys. Rev. D 42 (1990) 3613.
!  [5] L.K. Resvanis  (The NESTOR Collaboration),  "NESTOR: A Neutrino
!      Particle Astrophysics Underwater Laboratory for the  Miditerra-
!      nean," in Proc. of the Workshop  on High Energy Neutrino Astro-
!      physics, Honolulu, Hawaii, 1992,  edited  by V.J.Stenger et al.
!      (World Scientific, 1992), p.325.
!  [6] I.W. Rogers and M. Tristam, "The Absolute Depth-Intensity  Cur-
!      ve for Cosmic-Ray Muons  Underwater and the  Integral Sea-Level
!      Momentum Spectrum in the Range 1-100 GeV/c,"  J. Phys. G: Nucl.
!      Phys. 10 (1984) 983.
!  [7] I.A. Belolaptikov et al. (The BAIKAL Collaboration), "Status of
!      the Baikal Neutrino Detector," to appear in Proc. of the  26'th
!      International Conference on High Energy Physics, Dallas, 1992.
!  [8] C. Castagnoli and  O. Saavedra,  "Cosmic-Ray Muon Search at Mt.
!      Blanc Laboratory," Nuovo Cimento 9 C (1986) 111.
!  [9] H.E. Bergeson,  J.W. Keuffel, M.O. Larson, E.R. Martin and G.W.
!      Mason,  "Evidence for a new  production process  for $10^12$ eV
!      muons," Phys. Rev. Lett. 19 (1967) 1487.
! [10] M.F. Crouch, P.B. Landecker,  J.F. Lathrop et al.,  "Cosmic-ray
!      muon fluxes deep  underground: intensity vs depth and the neut-
!      rinoinduced component," Phys. Rev. D 18 (1978) 2239.       
! [11] S. Miyake,  V.S. Narasimham, and  P.V. Ramana Murthy,  Il Nuovo
!      Cimento 32 (1964) 1505;  M.R. Krishnaswamy,  M.G.K. Menon, V.S.
!      Narasimham et al., "Observations on the cosmic ray angular dis-
!      tribution underground at 1500 hg/cm$^2$ (Kolar) and the Questi-
!      on  of muon production in the TeV Energy Region," Phys.Lett. 27
!      B (1968) 535.
! [12] Ch. Berger, M. Frohlich, H. Monch et al. (The Frejus Collabora-
!      tion),  "Experimental  study of  muon bundles  observed  in the
!      Frejus detector," Phys. Rev. D 40 (1989) 2163.
! [13] S.P. Ahlen, M. Ambrosio, G. Auriemma et al. (The MACRO Collabo-
!      ration), "Study of penetrating cosmic ray  muons and search for
!      large scale anisotropies at the  Gran Sasso Laboratory,"  Phys.
!      Lett. 249B (1990) 149.
! [14] V.I. Gurentsov,  "Calculation of intensity end energy characte-
!      ristics of cosmic ray  muons in the place  of situation of  the
!      scintillation  telescope  BNO," INR P-0379  (Moscow, 1984)  [in
!      Russian]; Yu.M. Andreev, V.I. Gurentsov, and  I.M. Kogai, "Muon
!      Intensity from the Baksan Underground Scintillation Telescope,"
!      in Proc. of the 20'th International Cosmic Ray Conference, Mos-
!      cow, USSR, August 2 - 15, 1987, edited by  V.A. Kozyarivsky  et
!      al. (Nauka, Moscow, 1987), Vol.6, p.200.
! [15] P.T. Hayman  and A.W. Wolfendale,  Proc.  Roy. Soc. 275  (1963)
!      391.
!----------------------------------------------------------------------!
!    hnu_p=28.816*sqrt(rho*Z/A)                                        !Plasma energy [eV]
!        C=(2*log(I/hnu_p)+1)                                          !Generic formula
!      -------------------------- 
!     || Laboratory  |  I [eV]  ||
!     ||-------------|----------||
!     || Mont Blanc  |  141.0   ||
!     || KGF         |  165.0   || I=12.9*(9.76+58.8/(12.9**1.19))
!     || Baksan      |  155.3   ||
!     || Gran Sasso  |  151.8   ||
!     || Frejus      |  154.5   ||
!      -------------------------- 
!----------------------------------------------------------------------!
! Material constants for elemental media. Values are given for the ratio 
! of atomic  number-to-mass Z/A, the mean excitation energy I, and the 
! density.  Some density values are only nominal; those for Z = 85 and 
! 87 were arbitrarily set to 10 in order to complete the calculations. 
! The table is taken from 
! http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
!     -------------------------------------------------------------- 
!    |   Z      Element               Z/A      I[eV] Density[g/cm^3]|
!     -------------------------------------------------------------- 
!    |   1  H   HYDROGEN            0.99212     19.2    8.375E-05   |
!    |   2  He  HELIUM              0.49968     41.8    1.663E-04   |
!    |   3  Li  LITHIUM             0.43221     40.0    5.340E-01   |
!    |   4  Be  BERYLLIUM           0.44384     63.7    1.848E+00   |
!    |   5  B   BORON               0.46245     76.0    2.370E+00   |
!    |   6  C   CARBON, GRAPHITE    0.49954     78.0    1.700E+00   |
!    |   7  N   NITROGEN            0.49976     82.0    1.165E-03   |
!    |   8  O   OXYGEN              0.50002     95.0    1.332E-03   |
!    |   9  F   FLUORINE            0.47372    115.0    1.580E-03   |
!    |  10  Ne  NEON                0.49555    137.0    8.385E-04   |
!    |  11  Na  SODIUM              0.47847    149.0    9.710E-01   |
!    |  12  Mg  MAGNESIUM           0.49373    156.0    1.740E+00   |
!    |  13  Al  ALUMINUM            0.48181    166.0    2.699E+00   |
!    |  14  Si  SILICON             0.49848    173.0    2.330E+00   |
!    |  15  P   PHOSPHORUS          0.48428    173.0    2.200E+00   |
!    |  16  S   SULFUR              0.49897    180.0    2.000E+00   |
!    |  17  Cl  CHLORINE            0.47951    174.0    2.995E-03   |
!    |  18  Ar  ARGON               0.45059    188.0    1.662E-03   |
!    |  19  K   POTASSIUM           0.48595    190.0    8.620E-01   |
!    |  20  Ca  CALCIUM             0.49903    191.0    1.550E+00   |
!    |  21  Sc  SCANDIUM            0.46712    216.0    2.989E+00   |
!    |  22  Ti  TITANIUM            0.45948    233.0    4.540E+00   |
!    |  23  V   VANADIUM            0.45150    245.0    6.110E+00   |
!    |  24  Cr  CHROMIUM            0.46157    257.0    7.180E+00   |
!    |  25  Mn  MANGANESE           0.45506    272.0    7.440E+00   |
!    |  26  Fe  IRON                0.46556    286.0    7.874E+00   |
!    |  27  Co  COBALT              0.45815    297.0    8.900E+00   |
!    |  28  Ni  NICKEL              0.47708    311.0    8.902E+00   |
!    |  29  Cu  COPPER              0.45636    322.0    8.960E+00   |
!    |  30  Zn  ZINC                0.45879    330.0    7.133E+00   |
!    |  31  Ga  GALLIUM             0.44462    334.0    5.904E+00   |
!    |  32  Ge  GERMANIUM           0.44071    350.0    5.323E+00   |
!    |  33  As  ARSENIC             0.44046    347.0    5.730E+00   |
!    |  34  Se  SELENIUM            0.43060    348.0    4.500E+00   |
!    |  35  Br  BROMINE             0.43803    343.0    7.072E-03   |
!    |  36  Kr  KRYPTON             0.42959    352.0    3.478E-03   |
!    |  37  Rb  RUBIDIUM            0.43291    363.0    1.532E+00   |
!    |  38  Sr  STRONTIUM           0.43369    366.0    2.540E+00   |
!    |  39  Y   YTTRIUM             0.43867    379.0    4.469E+00   |
!    |  40  Zr  ZIRCONIUM           0.43848    393.0    6.506E+00   |
!    |  41  Nb  NIOBIUM             0.44130    417.0    8.570E+00   |
!    |  42  Mo  MOLYBDENUM          0.43777    424.0    1.022E+01   |
!    |  43  Tc  TECHNETIUM          0.43919    428.0    1.150E+01   |
!    |  44  Ru  RUTHENIUM           0.43534    441.0    1.241E+01   |
!    |  45  Rh  RHODIUM             0.43729    449.0    1.241E+01   |
!    |  46  Pd  PALLADIUM           0.43225    470.0    1.202E+01   |
!    |  47  Ag  SILVER              0.43572    470.0    1.050E+01   |
!    |  48  Cd  CADMIUM             0.42700    469.0    8.650E+00   |
!    |  49  In  INDIUM              0.42676    488.0    7.310E+00   |
!    |  50  Sn  TIN                 0.42120    488.0    7.310E+00   |
!    |  51  Sb  ANTIMONY            0.41889    487.0    6.691E+00   |
!    |  52  Te  TELLURIUM           0.40752    485.0    6.240E+00   |
!    |  53  I   IODINE              0.41764    491.0    4.930E+00   |
!    |  54  Xe  XENON               0.41130    482.0    5.485E-03   |
!    |  55  Cs  CESIUM              0.41383    488.0    1.873E+00   |
!    |  56  Ba  BARIUM              0.40779    491.0    3.500E+00   |
!    |  57  La  LANTHANUM           0.41035    501.0    6.154E+00   |
!    |  58  Ce  CERIUM              0.41395    523.0    6.657E+00   |
!    |  59  Pr  PRASEODYMIUM        0.41871    535.0    6.710E+00   |
!    |  60  Nd  NEODYMIUM           0.41597    546.0    6.900E+00   |
!    |  61  Pm  PROMETHIUM          0.42094    560.0    7.220E+00   |
!    |  62  Sm  SAMARIUM            0.41234    574.0    7.460E+00   |
!    |  63  Eu  EUROPIUM            0.41457    580.0    5.243E+00   |
!    |  64  Gd  GADOLINIUM          0.40699    591.0    7.900E+00   |
!    |  65  Tb  TERBIUM             0.40900    614.0    8.229E+00   |
!    |  66  Dy  DYSPROSIUM          0.40615    628.0    8.550E+00   |
!    |  67  Ho  HOLMIUM             0.40623    650.0    8.795E+00   |
!    |  68  Er  ERBIUM              0.40655    658.0    9.066E+00   |
!    |  69  Tm  THULIUM             0.40844    674.0    9.321E+00   |
!    |  70  Yb  YTTERBIUM           0.40453    684.0    6.730E+00   |
!    |  71  Lu  LUTETIUM            0.40579    694.0    9.840E+00   |
!    |  72  Hf  HAFNIUM             0.40338    705.0    1.331E+01   |
!    |  73  Ta  TANTALUM            0.40343    718.0    1.665E+01   |
!    |  74  W   TUNGSTEN            0.40250    727.0    1.930E+01   |
!    |  75  Re  RHENIUM             0.40278    736.0    2.102E+01   |
!    |  76  Os  OSMIUM              0.39958    746.0    2.257E+01   |
!    |  77  Ir  IRIDIUM             0.40058    757.0    2.242E+01   |
!    |  78  Pt  PLATINUM            0.39984    790.0    2.145E+01   |
!    |  79  Au  GOLD                0.40108    790.0    1.932E+01   |
!    |  80  Hg  MERCURY             0.39882    800.0    1.355E+01   |
!    |  81  Tl  THALLIUM            0.39631    810.0    1.172E+01   |
!    |  82  Pb  LEAD                0.39575    823.0    1.135E+01   |
!    |  83  Bi  BISMUTH             0.39717    823.0    9.747E+00   |
!    |  84  Po  POLONIUM            0.40195    830.0    9.320E+00   |
!    |  85  At  ASTATINE            0.40479    825.0    1.000E+01   |
!    |  86  Rn  RADON               0.38736    794.0    9.066E-03   |
!    |  87  Fr  FRANCIUM            0.39010    827.0    1.000E+01   |
!    |  88  Ra  RADIUM              0.38934    826.0    5.000E+00   |
!    |  89  Ac  ACTINIUM            0.39202    841.0    1.007E+01   |
!    |  90  Th  THORIUM             0.38787    847.0    1.172E+01   |
!    |  91  Pa  PROTACTINIUM        0.39388    878.0    1.537E+01   |
!    |  92  U   URANIUM             0.38651    890.0    1.895E+01   |
!     -------------------------------------------------------------- 
!**********************************************************************!
