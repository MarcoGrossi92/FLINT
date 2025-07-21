module ecker_mod
  implicit none
contains  
  subroutine ecker(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic
    use U_Lib_Chemistry_data
    implicit none
    real(8), intent(in)  :: roi(nsc)
    real(8), intent(in)  :: temp
    real(8), intent(out) :: omegadot(nsc) 
    real(8), intent(in)  :: rotot

    real(8) :: coi(nsc+1), Tdiff 
    real(8) :: M !< Third body
    integer :: is, T_i, Tint(2)
    real(8) :: prodf(1:28), prodb(1:28)

    do is = 1, nsc 
    coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
    enddo 
    T_i = int(temp) 
    Tdiff  = temp-T_i 
    Tint(1) = T_i 
    Tint(2) = T_i + 1 
    ! reac n. 1: H2 + O2 <=> H + HO2
    prodf(1)=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(2)**1.0)
    prodb(1)=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(3)**1.0)
    ! reac n. 2: H + O2 <=> O + OH
    prodf(2)=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(2)**1.0)
    prodb(2)=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(6)**1.0)
    ! reac n. 3: H2 + O <=> H + OH
    prodf(3)=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(5)**1.0)
    prodb(3)=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(6)**1.0)
    ! reac n. 4: H2 + OH <=> H + H2O
    prodf(4)=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(6)**1.0)
    prodb(4)=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(7)**1.0)
    ! reac n. 5: 2 OH <=> H2O + O
    prodf(5)=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*(coi(6)**2.0)
    prodb(5)=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(5)**1.0)
    ! reac n. 6: H + OH + M <=> H2O + M
    M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*6.0+coi(8)*1.0+coi(9)*1.0+coi(10)*1.0+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0
    prodf(6)=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(6)**1.0)*M
    prodb(6)=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*M
    ! reac n. 7: 2 H + M <=> H2 + M
    M=0.d0+coi(1)*2.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*6.0+coi(8)*1.0+coi(9)*1.0+coi(10)*1.0+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0
    prodf(7)=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*(coi(4)**2.0)*M
    prodb(7)=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*(coi(1)**1.0)*M
    ! reac n. 8: H + O + M <=> OH + M
    M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*5.0+coi(8)*1.0+coi(9)*1.0+coi(10)*1.0+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0
    prodf(8)=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(5)**1.0)*M
    prodb(8)=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*(coi(6)**1.0)*M
    ! reac n. 9: H + O2 + M <=> HO2 + M
    M=0.d0+coi(1)*2.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*16.0+coi(8)*1.0+coi(9)*1.0+coi(10)*1.0+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0
    prodf(9)=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(2)**1.0)*M
    prodb(9)=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*M
    ! reac n. 10: H + HO2 <=> 2 OH
    prodf(10)=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(3)**1.0)
    prodb(10)=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*(coi(6)**2.0)
    ! reac n. 11: H2 + HO2 <=> H2O + OH
    prodf(11)=comp_ch_tabT(11,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(3)**1.0)
    prodb(11)=comp_ch_tabT(11,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(6)**1.0)
    ! reac n. 12: HO2 + O <=> O2 + OH
    prodf(12)=comp_ch_tabT(12,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(5)**1.0)
    prodb(12)=comp_ch_tabT(12,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(6)**1.0)
    ! reac n. 13: HO2 + OH <=> H2O + O2
    prodf(13)=comp_ch_tabT(13,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(6)**1.0)
    prodb(13)=comp_ch_tabT(13,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(2)**1.0)
    ! reac n. 14: 2 HO2 <=> H2O2 + O2
    prodf(14)=comp_ch_tabT(14,kf_tab,Tint,Tdiff)*(coi(3)**2.0)
    prodb(14)=comp_ch_tabT(14,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(2)**1.0)
    ! reac n. 15: H + H2O2 <=> H2 + HO2
    prodf(15)=comp_ch_tabT(15,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(8)**1.0)
    prodb(15)=comp_ch_tabT(15,kb_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(3)**1.0)
    ! reac n. 16: H2O2 + O <=> HO2 + OH
    prodf(16)=comp_ch_tabT(16,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(5)**1.0)
    prodb(16)=comp_ch_tabT(16,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(6)**1.0)
    ! reac n. 17: H2O2 + OH <=> H2O + HO2
    prodf(17)=comp_ch_tabT(17,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(6)**1.0)
    prodb(17)=comp_ch_tabT(17,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(3)**1.0)
    ! reac n. 18: H2O2 + M <=> 2 OH + M
    M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*15.0+coi(8)*1.0+coi(9)*1.0+coi(10)*1.0+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0
    prodf(18)=comp_ch_tabT(18,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*M
    prodb(18)=comp_ch_tabT(18,kb_tab,Tint,Tdiff)*(coi(6)**2.0)*M
    ! reac n. 19: 2 O + M <=> O2 + M
    M=sum(coi(1:14))
    prodf(19)=comp_ch_tabT(19,kf_tab,Tint,Tdiff)*(coi(5)**2.0)*M
    prodb(19)=comp_ch_tabT(19,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*M
    ! reac n. 20: CO + OH <=> CO2 + H
    prodf(20)=comp_ch_tabT(20,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(6)**1.0)
    prodb(20)=comp_ch_tabT(20,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(4)**1.0)
    ! reac n. 21: CO + O2 <=> CO2 + O
    prodf(21)=comp_ch_tabT(21,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(2)**1.0)
    prodb(21)=comp_ch_tabT(21,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(5)**1.0)
    ! reac n. 22: CO + O + M <=> CO2 + M
    M=0.d0+coi(1)*2.5+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*12.0+coi(8)*1.0+coi(9)*1.9+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0
    prodf(22)=comp_ch_tabT(22,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(5)**1.0)*M
    prodb(22)=comp_ch_tabT(22,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*M
    ! reac n. 23: H + HCL <=> CL + H2
    prodf(23)=comp_ch_tabT(23,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(11)**1.0)
    prodb(23)=comp_ch_tabT(23,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(1)**1.0)
    ! reac n. 24: CL2 + H <=> CL + HCL
    prodf(24)=comp_ch_tabT(24,kf_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(4)**1.0)
    prodb(24)=comp_ch_tabT(24,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(11)**1.0)
    ! reac n. 25: HCL + OH <=> CL + H2O
    prodf(25)=comp_ch_tabT(25,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(6)**1.0)
    prodb(25)=comp_ch_tabT(25,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(7)**1.0)
    ! reac n. 26: HCL + O <=> CL + OH
    prodf(26)=comp_ch_tabT(26,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(5)**1.0)
    prodb(26)=comp_ch_tabT(26,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(6)**1.0)
    ! reac n. 27: 2 CL + M <=> CL2 + M
    M=sum(coi(1:14))
    prodf(27)=comp_ch_tabT(27,kf_tab,Tint,Tdiff)*(coi(13)**2.0)*M
    prodb(27)=comp_ch_tabT(27,kb_tab,Tint,Tdiff)*(coi(12)**1.0)*M
    ! reac n. 28: CL + H + M <=> HCL + M
    M=sum(coi(1:14))
    prodf(28)=comp_ch_tabT(28,kf_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(4)**1.0)*M
    prodb(28)=comp_ch_tabT(28,kb_tab,Tint,Tdiff)*(coi(11)**1.0)*M
    ! species source terms
    omegadot(1)=Wm_tab(1)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(23)-prodb(23)))
    omegadot(2)=Wm_tab(2)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(21)-prodb(21)))
    omegadot(3)=Wm_tab(3)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-2.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(17)-prodb(17)))
    omegadot(4)=Wm_tab(4)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(6)-prodb(6))+(0.0-2.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(28)-prodb(28)))
    omegadot(5)=Wm_tab(5)*(+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-2.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(26)-prodb(26)))
    omegadot(6)=Wm_tab(6)*(+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-2.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(8)-prodb(8))+(2.0-0.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(2.0-0.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26)))
    omegadot(7)=Wm_tab(7)*(+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(25)-prodb(25)))
    omegadot(8)=Wm_tab(8)*(+(1.0-0.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18)))
    omegadot(9)=Wm_tab(9)*(+(0.0-1.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22)))
    omegadot(10)=Wm_tab(10)*(+(1.0-0.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(22)-prodb(22)))
    omegadot(11)=Wm_tab(11)*(+(0.0-1.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(25)-prodb(25))+(0.0-1.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(28)-prodb(28)))
    omegadot(12)=Wm_tab(12)*(+(0.0-1.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(27)-prodb(27)))
    omegadot(13)=Wm_tab(13)*(+(1.0-0.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26))+(0.0-2.0)*(prodf(27)-prodb(27))+(0.0-1.0)*(prodf(28)-prodb(28)))
    omegadot(14)=0.d0
  end subroutine ecker
end module ecker_mod