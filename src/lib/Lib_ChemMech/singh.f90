module singh_mod
  implicit none
contains

  ! Singh (1994) global mechanism for ethylene combustion.
  ! 10 reversible reactions for 9 species.
  subroutine Singh(roi,temp,omegadot)
    use U_Lib_Thermodynamic
    use U_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(nsc)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(nsc)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: M, prod1f, prod1b, prod1, prod2f, prod2b, prod2
    real(8) :: prod3f, prod3b, prod3, prod4f, prod4b, prod4
    real(8) :: prod5f, prod5b, prod5, prod6f, prod6b, prod6
    real(8) :: prod7f, prod7b, prod7, prod8f, prod8b, prod8
    real(8) :: prod9f, prod9b, prod9, prod10f, prod10b, prod10

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 3rd body eff: 2.5 H2, 16.0 H2O, 1.0 else
    M = sum(coi(1:5))+2.5d0*coi(6)+16d0*coi(7)+coi(8)+coi(9)
    ! C2H4 + O2 <-> 2CO+ 2H2
    prod1f=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*coi(2)*coi(8)
    prod1b=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*coi(4)**2*coi(6)**2
    prod1=prod1f-prod1b
    ! CO + O <-> CO2 + M
    prod2f=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(4)*coi(9)*M
    prod2b=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*coi(5)*M
    prod2=prod2f-prod2b
    !CO + OH <-> CO2 + H
    prod3f=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(4)*coi(3)
    prod3b=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(5)*coi(1)
    prod3=prod3f-prod3b
    !H2 + O2 <-> OH+OH
    prod4f=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*coi(6)*coi(8)
    prod4b=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*coi(3)*coi(3)
    prod4=prod4f-prod4b
    !H + O2 <-> OH + O
    prod5f=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1)*coi(8)
    prod5b=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(3)*coi(9)
    prod5=prod5f-prod5b
    !OH +H2 <-> H2O + H
    prod6f=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3)*coi(6)
    prod6b=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(7)*coi(1)
    prod6=prod6f-prod6b
    !O + H2 <-> OH + H
    prod7f=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(9)*coi(6)
    prod7b=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(3)*coi(1)
    prod7=prod7f-prod7b
    !OH + OH <-> H2O + O
    prod8f=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*coi(3)*coi(3)
    prod8b=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*coi(7)*coi(9)
    prod8=prod8f-prod8b
    !H + H <-> H2 + M
    prod9f=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*coi(1)*coi(1)*M
    prod9b=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*coi(6)*M
    prod9=prod9f-prod9b
    !H + OH <-> H2O + M
    prod10f=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*coi(1)*coi(3)*M
    prod10b=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*coi(7)*M
    prod10=prod10f-prod10b

    ! H
    omegadot(1)=wm_tab(1)*(prod3+prod6+prod7-2*prod9-prod10)
    ! C2H4
    omegadot(2)=wm_tab(2)*(-prod1)
    ! OH
    omegadot(3)=wm_tab(3)*(-prod3+2*prod4+prod5-prod6+ &
                            prod7-2*prod8-prod10)
    ! CO
    omegadot(4)=wm_tab(4)*(2*prod1-prod2-prod3)
    ! CO2
    omegadot(5)=wm_tab(5)*(prod2+prod3)
    ! H2
    omegadot(6)=wm_tab(6)*(2*prod1-prod4-prod6-prod7+prod9)
    ! H2O
    omegadot(7)=wm_tab(7)*(prod6+prod8+prod10)
    ! O2
    omegadot(8)=wm_tab(8)*(-prod1-prod4-prod5)
    ! O
    omegadot(9)=wm_tab(9)*(-prod2+prod5-prod7+prod8)

  end subroutine Singh


  ! Singh (1994) + paraffin cracking
  subroutine Singh_WC32(roi,temp,omegadot)
    use U_Lib_Thermodynamic
    use U_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(nsc)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(nsc)
    ! Local
    integer :: is, T_i, Tint(2)
    integer, parameter :: iH=1, iC2H4=2, iOH=3, iCO=4, iCO2=5
    integer, parameter :: iH2=6, iH2O=7, iO2=8, iO=9, iC32H66=10
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: M, prodf(11), prodb(11), prod(11)

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 3rd body eff: 2.5 H2, 16.0 H2O, 1.0 else
    M = sum(coi(iH:iCO2)) + 2.5d0*coi(iH2) + 16d0*coi(iH2O) + coi(iO2) + coi(iO)
    ! C2H4 + O2 <-> 2CO+ 2H2
    prodf(1)=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*coi(iC2H4)*coi(iO2)
    prodb(1)=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*coi(iCO)**2*coi(iH2)**2
    prod(1)=prodf(1)-prodb(1)
    ! CO + O + M <-> CO2 + M
    prodf(2)=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iO)*M
    prodb(2)=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*coi(iCO2)*M
    prod(2)=prodf(2)-prodb(2)
    !CO + OH <-> CO2 + H
    prodf(3)=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iOH)
    prodb(3)=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(iCO2)*coi(iH)
    prod(3)=prodf(3)-prodb(3)
    !H2 + O2 <-> OH+OH
    prodf(4)=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*coi(iH2)*coi(iO2)
    prodb(4)=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*coi(iOH)**2
    prod(4)=prodf(4)-prodb(4)
    !H + O2 <-> OH + O
    prodf(5)=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(iH)*coi(iO2)
    prodb(5)=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iO)
    prod(5)=prodf(5)-prodb(5)
    !OH +H2 <-> H2O + H
    prodf(6)=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(iOH)*coi(iH2)
    prodb(6)=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iH)
    prod(6)=prodf(6)-prodb(6)
    !O + H2 <-> OH + H
    prodf(7)=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(iO)*coi(iH2)
    prodb(7)=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iH)
    prod(7)=prodf(7)-prodb(7)
    !OH + OH <-> H2O + O
    prodf(8)=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*coi(iOH)**2
    prodb(8)=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iO)
    prod(8)=prodf(8)-prodb(8)
    !H + H + M <-> H2 + M
    prodf(9)=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*coi(iH)**2*M
    prodb(9)=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*coi(iH2)*M
    prod(9)=prodf(9)-prodb(9)
    !H + OH <-> H2O + M
    prodf(10)=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*coi(iH)*coi(iOH)*M
    prodb(10)=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*coi(iH2O)*M
    prod(10)=prodf(10)-prodb(10)
    !C32H66 --> 16C2H4 + H2
    if (coi(iC32H66)>1d-10) then
      prodf(11)=comp_ch_tabT(11,kf_tab,Tint,Tdiff)*(coi(iC32H66))
    else
      prodf(11)=0d0
    endif
    prodb(11)=0d0
    prod(11)=prodf(11)

    ! H
    omegadot(iH)=wm_tab(iH)*(prod(3)-prod(5)+prod(6)+prod(7)-2*prod(9)-prod(10))
    ! C2H4
    omegadot(iC2H4)=wm_tab(iC2H4)*(-prod(1)+16*prod(11))
    ! OH
    omegadot(iOH)=wm_tab(iOH)*(-prod(3)+2*prod(4)+prod(5)-prod(6)+ &
                                prod(7)-2*prod(8)-prod(10))
    ! CO
    omegadot(iCO)=wm_tab(iCO)*(2*prod(1)-prod(2)-prod(3))
    ! CO2
    omegadot(iCO2)=wm_tab(iCO2)*(prod(2)+prod(3))
    ! H2
    omegadot(iH2)=wm_tab(iH2)*(2*prod(1)-prod(4)-prod(6)-prod(7)+prod(9)+prod(11))
    ! H2O
    omegadot(iH2O)=wm_tab(iH2O)*(prod(6)+prod(8)+prod(10))
    ! O2
    omegadot(iO2)=wm_tab(iO2)*(-prod(1)-prod(4)-prod(5))
    ! O
    omegadot(iO)=wm_tab(iO)*(-prod(2)+prod(5)-prod(7)+prod(8))
    ! C32H66
    omegadot(iC32H66)=wm_tab(iC32H66)*(-prod(11))

  end subroutine Singh_WC32

end module singh_mod