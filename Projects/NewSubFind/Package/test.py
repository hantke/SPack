import SPack_NewSubfind_iso_Polar as SP_NSF_iso
import SPack_NewSubfind_Polar as SP_NSF
f = 0.05
M_c = 7e11
e = 0.1
a = 0
djmin = 0.0
print SP_NSF.Tdf(M_c, M_c*0.3, djmin,e,a, 1e-4, 16, 0)
#print SP_NSF_iso.Tdf(M_c, M_c*0.3, djmin,e,a, 1e-4, 16, 0)
#print SP_NSF.Tdf(M_c, M_c*0.3, djmin,e,a, 1e-5, 16, 1)/SP_NSF.Tdf(M_c, M_c*f, djmin,e,a, 1e-5, 16, 1)
