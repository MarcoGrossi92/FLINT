module Gerlinger9_mod
implicit none
contains
subroutine Gerlinger9(roi,temp,omegadot)
use FLINT_Lib_Thermodynamic
use FLINT_Lib_Chemistry_data
use FLINT_Lib_Chemistry_falloff
implicit none
real(8), intent(inout)  :: roi(ns)
real(8), intent(in)  :: temp
real(8), intent(out) :: omegadot(ns) 

real(8) :: coi(ns), Tdiff 
real(8) :: M !< Third body
integer :: is, T_i, Tint(2)
real(8) :: prodf(1:19), prodb(1:19)
real(8) :: k(2) !< Falloff rate coefficients


do is = 1, ns 
 coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
enddo 
T_i = int(temp) 
Tdiff  = temp-T_i 
Tint(1) = T_i 
Tint(2) = T_i + 1 
! reac n. 1: H2 + O2 <=> H + HO2
prodf(1)=f_kf(1,Tint,Tdiff)*(coi(3))*(coi(2))
prodb(1)=f_kb(1,Tint,Tdiff)*(coi(7))*(coi(8))
! reac n. 2: H + O2 <=> O + OH
prodf(2)=f_kf(2,Tint,Tdiff)*(coi(7))*(coi(2))
prodb(2)=f_kb(2,Tint,Tdiff)*(coi(6))*(coi(5))
! reac n. 3: H2 + O <=> H + OH
prodf(3)=f_kf(3,Tint,Tdiff)*(coi(3))*(coi(6))
prodb(3)=f_kb(3,Tint,Tdiff)*(coi(7))*(coi(5))
! reac n. 4: H2 + OH <=> H + H2O
prodf(4)=f_kf(4,Tint,Tdiff)*(coi(3))*(coi(5))
prodb(4)=f_kb(4,Tint,Tdiff)*(coi(7))*(coi(4))
! reac n. 5: 2 OH <=> H2O + O
prodf(5)=f_kf(5,Tint,Tdiff)*(coi(5)*coi(5))
prodb(5)=f_kb(5,Tint,Tdiff)*(coi(4))*(coi(6))
! reac n. 6: H + OH + M <=> H2O + M
M=coi(1)+coi(2)+coi(3)+coi(4)*6+coi(5)+coi(6)+coi(7)+coi(8)+coi(9)
prodf(6)=f_kf(6,Tint,Tdiff)*(coi(7))*(coi(5))*M
prodb(6)=f_kb(6,Tint,Tdiff)*(coi(4))*M
! reac n. 7: 2 H + M <=> H2 + M
M=coi(1)+coi(2)+coi(3)*2+coi(4)*6+coi(5)+coi(6)+coi(7)+coi(8)+coi(9)
prodf(7)=f_kf(7,Tint,Tdiff)*(coi(7)*coi(7))*M
prodb(7)=f_kb(7,Tint,Tdiff)*(coi(3))*M
! reac n. 8: H + O + M <=> OH + M
M=coi(1)+coi(2)+coi(3)+coi(4)*5+coi(5)+coi(6)+coi(7)+coi(8)+coi(9)
prodf(8)=f_kf(8,Tint,Tdiff)*(coi(7))*(coi(6))*M
prodb(8)=f_kb(8,Tint,Tdiff)*(coi(5))*M
! reac n. 9: H + O2 + M <=> HO2 + M
M=coi(1)+coi(2)+coi(3)*2+coi(4)*16+coi(5)+coi(6)+coi(7)+coi(8)+coi(9)
prodf(9)=f_kf(9,Tint,Tdiff)*(coi(7))*(coi(2))*M
prodb(9)=f_kb(9,Tint,Tdiff)*(coi(8))*M
! reac n. 10: H + HO2 <=> 2 OH
prodf(10)=f_kf(10,Tint,Tdiff)*(coi(7))*(coi(8))
prodb(10)=f_kb(10,Tint,Tdiff)*(coi(5)*coi(5))
! reac n. 11: H + HO2 <=> H2O + O
prodf(11)=f_kf(11,Tint,Tdiff)*(coi(7))*(coi(8))
prodb(11)=f_kb(11,Tint,Tdiff)*(coi(4))*(coi(6))
! reac n. 12: HO2 + O <=> O2 + OH
prodf(12)=f_kf(12,Tint,Tdiff)*(coi(8))*(coi(6))
prodb(12)=f_kb(12,Tint,Tdiff)*(coi(2))*(coi(5))
! reac n. 13: HO2 + OH <=> H2O + O2
prodf(13)=f_kf(13,Tint,Tdiff)*(coi(8))*(coi(5))
prodb(13)=f_kb(13,Tint,Tdiff)*(coi(4))*(coi(2))
! reac n. 14: 2 HO2 <=> H2O2 + O2
prodf(14)=f_kf(14,Tint,Tdiff)*(coi(8)*coi(8))
prodb(14)=f_kb(14,Tint,Tdiff)*(coi(9))*(coi(2))
! reac n. 15: H + H2O2 <=> H2 + HO2
prodf(15)=f_kf(15,Tint,Tdiff)*(coi(7))*(coi(9))
prodb(15)=f_kb(15,Tint,Tdiff)*(coi(3))*(coi(8))
! reac n. 16: H2O2 + O <=> HO2 + OH
prodf(16)=f_kf(16,Tint,Tdiff)*(coi(9))*(coi(6))
prodb(16)=f_kb(16,Tint,Tdiff)*(coi(8))*(coi(5))
! reac n. 17: H2O2 + OH <=> H2O + HO2
prodf(17)=f_kf(17,Tint,Tdiff)*(coi(9))*(coi(5))
prodb(17)=f_kb(17,Tint,Tdiff)*(coi(4))*(coi(8))
! reac n. 18: H2O2 + M <=> 2 OH + M
M=coi(1)+coi(2)+coi(3)+coi(4)*15+coi(5)+coi(6)+coi(7)+coi(8)+coi(9)
prodf(18)=f_kf(18,Tint,Tdiff)*(coi(9))*M
prodb(18)=f_kb(18,Tint,Tdiff)*(coi(5)*coi(5))*M
! reac n. 19: 2 O + M <=> O2 + M
M=sum(coi(1:9))
prodf(19)=f_kf(19,Tint,Tdiff)*(coi(6)*coi(6))*M
prodb(19)=f_kb(19,Tint,Tdiff)*(coi(2))*M
! species source terms
omegadot(1)=0.d0
omegadot(2)=Wm_tab(2)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(19)-prodb(19)))
omegadot(3)=Wm_tab(3)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(15)-prodb(15)))
omegadot(4)=Wm_tab(4)*(+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(17)-prodb(17)))
omegadot(5)=Wm_tab(5)*(+(1.0-0.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-2.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(8)-prodb(8))+(2.0-0.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(2.0-0.0)*(prodf(18)-prodb(18)))
omegadot(6)=Wm_tab(6)*(+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-2.0)*(prodf(19)-prodb(19)))
omegadot(7)=Wm_tab(7)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(6)-prodb(6))+(0.0-2.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(15)-prodb(15)))
omegadot(8)=Wm_tab(8)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(13)-prodb(13))+(0.0-2.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(17)-prodb(17)))
omegadot(9)=Wm_tab(9)*(+(1.0-0.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18)))
end subroutine Gerlinger9
end module Gerlinger9_mod
