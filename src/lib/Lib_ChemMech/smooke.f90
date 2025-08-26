module smooke_mod
  implicit none
contains  
  subroutine smooke(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout)  :: roi(nsc)
    real(8), intent(in)  :: temp
    real(8), intent(out) :: omegadot(nsc) 

    real(8) :: coi(nsc+1), Tdiff 
    real(8) :: M !< Third body
    integer :: is, T_i, Tint(2)
    real(8) :: prodf(1:35), prodb(1:35)

    do is = 1, nsc
    coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
    enddo 
    T_i = int(temp) 
    Tdiff  = temp-T_i 
    Tint(1) = T_i 
    Tint(2) = T_i + 1 
    ! reac n. 1: H + O2 => O + OH
    prodf(1)=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(3)**1.0)
    prodb(1)=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(6)**1.0)
    ! reac n. 2: O + OH => H + O2
    prodf(2)=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(6)**1.0)
    prodb(2)=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(3)**1.0)
    ! reac n. 3: H2 + O => H + OH
    prodf(3)=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(4)**1.0)
    prodb(3)=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(6)**1.0)
    ! reac n. 4: H + OH => H2 + O
    prodf(4)=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(6)**1.0)
    prodb(4)=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(4)**1.0)
    ! reac n. 5: H2 + OH => H + H2O
    prodf(5)=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(6)**1.0)
    prodb(5)=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(9)**1.0)
    ! reac n. 6: H + H2O => H2 + OH
    prodf(6)=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(9)**1.0)
    prodb(6)=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(6)**1.0)
    ! reac n. 7: 2 OH => H2O + O
    prodf(7)=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*(coi(6)**2.0)
    prodb(7)=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(4)**1.0)
    ! reac n. 8: H2O + O => 2 OH
    prodf(8)=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(4)**1.0)
    prodb(8)=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*(coi(6)**2.0)
    ! reac n. 9: H + O2 + M => HO2 + M
    M=0.d0+coi(1)*6.5+coi(2)*1.0+coi(3)*0.4+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*1.0+coi(8)*1.0+coi(9)*6.5+coi(10)*0.75+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.5+coi(16)*0.4
    prodf(9)=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(3)**1.0)*M
    prodb(9)=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*M
    ! reac n. 10: H + HO2 => 2 OH
    prodf(10)=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(7)**1.0)
    prodb(10)=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*(coi(6)**2.0)
    ! reac n. 11: H + HO2 => H2 + O2
    prodf(11)=comp_ch_tabT(11,kf_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(7)**1.0)
    prodb(11)=comp_ch_tabT(11,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(3)**1.0)
    ! reac n. 12: HO2 + OH => H2O + O2
    prodf(12)=comp_ch_tabT(12,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(6)**1.0)
    prodb(12)=comp_ch_tabT(12,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(3)**1.0)
    ! reac n. 13: CO + OH => CO2 + H
    prodf(13)=comp_ch_tabT(13,kf_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(6)**1.0)
    prodb(13)=comp_ch_tabT(13,kb_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(5)**1.0)
    ! reac n. 14: CO2 + H => CO + OH
    prodf(14)=comp_ch_tabT(14,kf_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(5)**1.0)
    prodb(14)=comp_ch_tabT(14,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(6)**1.0)
    ! reac n. 15: CH4 => CH3 + H
    prodf(15)=comp_ch_tabT(15,kf_tab,Tint,Tdiff)*(coi(1)**1.0)
    prodb(15)=comp_ch_tabT(15,kb_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(5)**1.0)
    ! reac n. 16: CH3 + H => CH4
    prodf(16)=comp_ch_tabT(16,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(5)**1.0)
    prodb(16)=comp_ch_tabT(16,kb_tab,Tint,Tdiff)*(coi(1)**1.0)
    ! reac n. 17: CH4 + H => CH3 + H2
    prodf(17)=comp_ch_tabT(17,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(5)**1.0)
    prodb(17)=comp_ch_tabT(17,kb_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(2)**1.0)
    ! reac n. 18: CH3 + H2 => CH4 + H
    prodf(18)=comp_ch_tabT(18,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(2)**1.0)
    prodb(18)=comp_ch_tabT(18,kb_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(5)**1.0)
    ! reac n. 19: CH4 + OH => CH3 + H2O
    prodf(19)=comp_ch_tabT(19,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(6)**1.0)
    prodb(19)=comp_ch_tabT(19,kb_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(9)**1.0)
    ! reac n. 20: CH3 + H2O => CH4 + OH
    prodf(20)=comp_ch_tabT(20,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(9)**1.0)
    prodb(20)=comp_ch_tabT(20,kb_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(6)**1.0)
    ! reac n. 21: CH3 + O => CH2O + H
    prodf(21)=comp_ch_tabT(21,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(4)**1.0)
    prodb(21)=comp_ch_tabT(21,kb_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(5)**1.0)
    ! reac n. 22: CH2O + H => H2 + HCO
    prodf(22)=comp_ch_tabT(22,kf_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(5)**1.0)
    prodb(22)=comp_ch_tabT(22,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(13)**1.0)
    ! reac n. 23: CH2O + OH => H2O + HCO
    prodf(23)=comp_ch_tabT(23,kf_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(6)**1.0)
    prodb(23)=comp_ch_tabT(23,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(13)**1.0)
    ! reac n. 24: H + HCO => CO + H2
    prodf(24)=comp_ch_tabT(24,kf_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(13)**1.0)
    prodb(24)=comp_ch_tabT(24,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(2)**1.0)
    ! reac n. 25: HCO + M => CO + H + M
    M=0.d0+coi(1)*6.5+coi(2)*1.0+coi(3)*0.4+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*1.0+coi(8)*1.0+coi(9)*6.5+coi(10)*0.75+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.5+coi(16)*0.4
    prodf(25)=comp_ch_tabT(25,kf_tab,Tint,Tdiff)*(coi(13)**1.0)*M
    prodb(25)=comp_ch_tabT(25,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(5)**1.0)*M
    ! reac n. 26: CH3 + O2 => CH3O + O
    prodf(26)=comp_ch_tabT(26,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(3)**1.0)
    prodb(26)=comp_ch_tabT(26,kb_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(4)**1.0)
    ! reac n. 27: CH3O + H => CH2O + H2
    prodf(27)=comp_ch_tabT(27,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(5)**1.0)
    prodb(27)=comp_ch_tabT(27,kb_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(2)**1.0)
    ! reac n. 28: CH3O + M => CH2O + H + M
    M=0.d0+coi(1)*6.5+coi(2)*1.0+coi(3)*0.4+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*1.0+coi(8)*1.0+coi(9)*6.5+coi(10)*0.75+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.5+coi(16)*0.4
    prodf(28)=comp_ch_tabT(28,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*M
    prodb(28)=comp_ch_tabT(28,kb_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(5)**1.0)*M
    ! reac n. 29: 2 HO2 => H2O2 + O2
    prodf(29)=comp_ch_tabT(29,kf_tab,Tint,Tdiff)*(coi(7)**2.0)
    prodb(29)=comp_ch_tabT(29,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(3)**1.0)
    ! reac n. 30: H2O2 + M => 2 OH + M
    M=0.d0+coi(1)*6.5+coi(2)*1.0+coi(3)*0.4+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*1.0+coi(8)*1.0+coi(9)*6.5+coi(10)*0.75+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.5+coi(16)*0.4
    prodf(30)=comp_ch_tabT(30,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*M
    prodb(30)=comp_ch_tabT(30,kb_tab,Tint,Tdiff)*(coi(6)**2.0)*M
    ! reac n. 31: 2 OH + M => H2O2 + M
    M=0.d0+coi(1)*6.5+coi(2)*1.0+coi(3)*0.4+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*1.0+coi(8)*1.0+coi(9)*6.5+coi(10)*0.75+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.5+coi(16)*0.4
    prodf(31)=comp_ch_tabT(31,kf_tab,Tint,Tdiff)*(coi(6)**2.0)*M
    prodb(31)=comp_ch_tabT(31,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*M
    ! reac n. 32: H2O2 + OH => H2O + HO2
    prodf(32)=comp_ch_tabT(32,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(6)**1.0)
    prodb(32)=comp_ch_tabT(32,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(7)**1.0)
    ! reac n. 33: H2O + HO2 => H2O2 + OH
    prodf(33)=comp_ch_tabT(33,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(7)**1.0)
    prodb(33)=comp_ch_tabT(33,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(6)**1.0)
    ! reac n. 34: H + OH + M => H2O + M
    M=0.d0+coi(1)*6.5+coi(2)*1.0+coi(3)*0.4+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*1.0+coi(8)*1.0+coi(9)*6.5+coi(10)*0.75+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.5+coi(16)*0.4
    prodf(34)=comp_ch_tabT(34,kf_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(6)**1.0)*M
    prodb(34)=comp_ch_tabT(34,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*M
    ! reac n. 35: 2 H + M => H2 + M
    M=0.d0+coi(1)*6.5+coi(2)*1.0+coi(3)*0.4+coi(4)*1.0+coi(5)*1.0+coi(6)*1.0+coi(7)*1.0+coi(8)*1.0+coi(9)*6.5+coi(10)*0.75+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.5+coi(16)*0.4
    prodf(35)=comp_ch_tabT(35,kf_tab,Tint,Tdiff)*(coi(5)**2.0)*M
    prodb(35)=comp_ch_tabT(35,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*M
    ! species source terms
    omegadot(1)=Wm_tab(1)*(+(0.0-1.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(20)-prodb(20)))
    omegadot(2)=Wm_tab(2)*(+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(35)-prodb(35)))
    omegadot(3)=Wm_tab(3)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(29)-prodb(29)))
    omegadot(4)=Wm_tab(4)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(26)-prodb(26)))
    omegadot(5)=Wm_tab(5)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(25)-prodb(25))+(0.0-1.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(34)-prodb(34))+(0.0-2.0)*(prodf(35)-prodb(35)))
    omegadot(6)=Wm_tab(6)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(0.0-2.0)*(prodf(7)-prodb(7))+(2.0-0.0)*(prodf(8)-prodb(8))+(2.0-0.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(19)-prodb(19))+(1.0-0.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(23)-prodb(23))+(2.0-0.0)*(prodf(30)-prodb(30))+(0.0-2.0)*(prodf(31)-prodb(31))+(0.0-1.0)*(prodf(32)-prodb(32))+(1.0-0.0)*(prodf(33)-prodb(33))+(0.0-1.0)*(prodf(34)-prodb(34)))
    omegadot(7)=Wm_tab(7)*(+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-2.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(32)-prodb(32))+(0.0-1.0)*(prodf(33)-prodb(33)))
    omegadot(8)=Wm_tab(8)*(+(1.0-0.0)*(prodf(29)-prodb(29))+(0.0-1.0)*(prodf(30)-prodb(30))+(1.0-0.0)*(prodf(31)-prodb(31))+(0.0-1.0)*(prodf(32)-prodb(32))+(1.0-0.0)*(prodf(33)-prodb(33)))
    omegadot(9)=Wm_tab(9)*(+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(32)-prodb(32))+(0.0-1.0)*(prodf(33)-prodb(33))+(1.0-0.0)*(prodf(34)-prodb(34)))
    omegadot(10)=Wm_tab(10)*(+(0.0-1.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(25)-prodb(25)))
    omegadot(11)=Wm_tab(11)*(+(1.0-0.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(26)-prodb(26)))
    omegadot(12)=Wm_tab(12)*(+(1.0-0.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(28)-prodb(28)))
    omegadot(13)=Wm_tab(13)*(+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(25)-prodb(25)))
    omegadot(14)=Wm_tab(14)*(+(1.0-0.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(27)-prodb(27))+(0.0-1.0)*(prodf(28)-prodb(28)))
    omegadot(15)=Wm_tab(15)*(+(1.0-0.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14)))
    omegadot(16)=0.d0
  end subroutine smooke
end module smooke_mod