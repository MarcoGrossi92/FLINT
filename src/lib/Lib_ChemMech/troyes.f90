module troyes_mod
  implicit none
contains
    subroutine troyes(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout)  :: roi(ns)
    real(8), intent(in)  :: temp
    real(8), intent(out) :: omegadot(ns) 

    real(8) :: coi(ns), Tdiff 
    real(8) :: M !< Third body
    integer :: is, T_i, Tint(2)
    real(8) :: prodf(1:17), prodb(1:17)

    do is = 1, ns
    roi(is) = max(roi(is), 0.d0)
    coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
    enddo 
    T_i = int(temp) 
    Tdiff  = temp-T_i 
    Tint(1) = T_i 
    Tint(2) = T_i + 1 
    ! reac n. 1: H + O2 <=> O + OH
    prodf(1)=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(2)**1.0)
    prodb(1)=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(3)**1.0)
    ! reac n. 2: H2 + O <=> H + OH
    prodf(2)=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(4)**1.0)
    prodb(2)=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(3)**1.0)
    ! reac n. 3: H2 + OH <=> H + H2O
    prodf(3)=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(3)**1.0)
    prodb(3)=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(5)**1.0)
    ! reac n. 4: 2 OH <=> H2O + O
    prodf(4)=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(3)**2.0)
    prodb(4)=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(4)**1.0)
    ! reac n. 5: 2 H + M <=> H2 + M
    M=sum(coi(1:12))
    prodf(5)=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*(coi(1)**2.0)*M
    prodb(5)=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*(coi(6)**1.0)*M
    ! reac n. 6: H + OH + M <=> H2O + M
    M=sum(coi(1:12))
    prodf(6)=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(3)**1.0)*M
    prodb(6)=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*(coi(5)**1.0)*M
    ! reac n. 7: H + O + M <=> OH + M
    M=sum(coi(1:12))
    prodf(7)=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(4)**1.0)*M
    prodb(7)=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*M
    ! reac n. 8: 2 O + M <=> O2 + M
    M=sum(coi(1:12))
    prodf(8)=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*(coi(4)**2.0)*M
    prodb(8)=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*M
    ! reac n. 9: CO + OH <=> CO2 + H
    prodf(9)=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(3)**1.0)
    prodb(9)=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(1)**1.0)
    ! reac n. 10: CO + O2 <=> CO2 + O
    prodf(10)=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(2)**1.0)
    prodb(10)=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(4)**1.0)
    ! reac n. 11: CO + O + M <=> CO2 + M
    M=sum(coi(1:12))
    prodf(11)=comp_ch_tabT(11,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(4)**1.0)*M
    prodb(11)=comp_ch_tabT(11,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*M
    ! reac n. 12: H + HCL <=> CL + H2
    prodf(12)=comp_ch_tabT(12,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*(coi(11)**1.0)
    prodb(12)=comp_ch_tabT(12,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(6)**1.0)
    ! reac n. 13: CL2 + H <=> CL + HCL
    prodf(13)=comp_ch_tabT(13,kf_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(1)**1.0)
    prodb(13)=comp_ch_tabT(13,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(11)**1.0)
    ! reac n. 14: HCL + OH <=> CL + H2O
    prodf(14)=comp_ch_tabT(14,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(3)**1.0)
    prodb(14)=comp_ch_tabT(14,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(5)**1.0)
    ! reac n. 15: HCL + O <=> CL + OH
    prodf(15)=comp_ch_tabT(15,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(4)**1.0)
    prodb(15)=comp_ch_tabT(15,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(3)**1.0)
    ! reac n. 16: 2 CL + M <=> CL2 + M
    M=sum(coi(1:12))
    prodf(16)=comp_ch_tabT(16,kf_tab,Tint,Tdiff)*(coi(9)**2.0)*M
    prodb(16)=comp_ch_tabT(16,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*M
    ! reac n. 17: CL + H + M <=> HCL + M
    M=sum(coi(1:12))
    prodf(17)=comp_ch_tabT(17,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(1)**1.0)*M
    prodb(17)=comp_ch_tabT(17,kb_tab,Tint,Tdiff)*(coi(11)**1.0)*M
    ! species source terms
    omegadot(1)=Wm_tab(1)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-2.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(17)-prodb(17)))
    omegadot(2)=Wm_tab(2)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(10)-prodb(10)))
    omegadot(3)=Wm_tab(3)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(0.0-2.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(15)-prodb(15)))
    omegadot(4)=Wm_tab(4)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(7)-prodb(7))+(0.0-2.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(15)-prodb(15)))
    omegadot(5)=Wm_tab(5)*(+(1.0-0.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(14)-prodb(14)))
    omegadot(6)=Wm_tab(6)*(+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(12)-prodb(12)))
    omegadot(7)=Wm_tab(7)*(+(0.0-1.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11)))
    omegadot(8)=Wm_tab(8)*(+(1.0-0.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11)))
    omegadot(9)=Wm_tab(9)*(+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(15)-prodb(15))+(0.0-2.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17)))
    omegadot(10)=Wm_tab(10)*(+(0.0-1.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(16)-prodb(16)))
    omegadot(11)=Wm_tab(11)*(+(0.0-1.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(17)-prodb(17)))
    omegadot(12)=0.d0
    end subroutine troyes
end module troyes_mod
