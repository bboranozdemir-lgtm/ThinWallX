# ThinWallX Phase v0.7 — Bileşik Gerilme Analizi, Bimoment Mekaniği ve Tepe Noktası Taraması (Stress Recovery & Non-Uniform Torsion)

## 1. Giriş ve Kapsam
v0.1–v0.6 fazlarında dondurulan rijitlik karakterizasyonunu (A, I, J, S, omega*, Cw), kesite uygulanan dış yükler altında (N, Vx, Vy, Mx, My, Tsv, B, M_omega) sac liflerinde analitik gerilme geri kazanımına (stress recovery) dönüştürmektir.
FEM veya 2D diskretizasyon kesinlikle kullanılmaz; tüm gerilme dağılımları orta-hat boyunca 1D analitik kapalı formüller ve polinom kökleriyle hesaplanır.

## 2. Matematiksel Temeller ve İşaret Sözleşmesi
Referans: Ağırlık merkezi C=(xc, yc), genel eksenlerde (Ix, Iy, Ixy). Kutup: Kayma merkezi S=(xs, ys).

### 2.1 Eksenel Normal Gerilme (sigma_zz(s))
Dış yükler: Eksenel kuvvet N [F], eğilme momentleri Mx, My [F*L] (sağ el kuralı) ve Vlasov Bimomenti B [F*L^2].
delta_I = Ix * Iy - Ixy^2 > 0 olmak üzere:
sigma_zz(s) = (N / A) 
            + ((My * Ix + Mx * Ixy) / delta_I) * (x(s) - xc) 
            - ((Mx * Iy + My * Ixy) / delta_I) * (y(s) - yc) 
            + (B * omega*(s)) / Cw

Özel Kurallar:
- Cw = 0 olan kesitlerde B != 0 uygulanırsa derhal `SingularSectionError("Section has zero warping resistance (Cw=0). Bimoment B cannot be sustained.")` fırlatılır. B == 0 ise terim tam 0.0'dır.
- sigma_zz(s) segment boyunca 1. dereceden bir polinomdur (doğrusaldır).

### 2.2 Kayma Akışı ve Kayma Gerilmesi Bileşenleri
Dış yükler: Enine kesme kuvvetleri Vx, Vy [F], Saint-Venant burulması Tsv [F*L] ve Çarpılma burulma momenti M_omega [F*L].

1. Enine Kesme Kayması (tau_V(s)):
   tau_V(s) = q_V(s) / t(s)  (v0.2/v0.5/v0.6 çözümleri, segment boyunca 2. derece parabolik).

2. Saint-Venant Serbest Burulma Kayması (tau_sv):
   Ortak burulma hızı theta' = Tsv / (G * J_total) kinematik sürekliliği uyarınca payda DAİMA J_total'dir:
   - Kapalı hücre kenarlarında (closed_edges): Membran akışı mevcuttur:
     tau_sv_membrane(s) = (Tsv / J_total) * (F_e / t_e)
     tau_sv_surface(s) = 0.0
   - Açık kollarda (open_edges): Membran akısı sıfırdır, sac yüzeyinde burulma kayması oluşur:
     tau_sv_membrane(s) = 0.0
     tau_sv_surface(s) = (abs(Tsv) / J_total) * t_e

3. İkincil Çarpılma Kayması (tau_omega(s)):
   Birim: M_omega birimi [F*L] (moment) olup akı birimi [F/L] üretir:
   Sektöriyel Statik Moment: S_omega(s) = integral_0^s omega*(xi) * t(xi) dxi
   Açık kollarda serbest uçta S_omega = 0.
   Kapalı hücrelerde oint (q_omega / t) ds = 0 uyumluluk döngüsü çözülerek q0_omega bulunur.
   q_omega(s) = - (M_omega / Cw) * S_omega_compatible(s)
   tau_omega(s) = q_omega(s) / t(s)  (segment boyunca 2. derece parabolik).

Net Süperpozisyonlar:
- Membran Kayması: tau_membrane(s) = tau_V(s) + tau_omega(s) + tau_sv_membrane(s)
- Kritik Dış Yüzey Kayması (Peak Envelope):
  tau_surface(s) = abs(tau_membrane(s)) + tau_sv_surface(s)

### 2.3 Von Mises Eşdeğer Gerilmesi, Ekstremum Taraması ve Kritik Yük Faktörü
Sac yüzeyindeki kritik eşdeğer gerilme:
sigma_vm(s) = sqrt( sigma_zz(s)^2 + 3 * (tau_surface(s))^2 )

Ekstremum Belirleme:
- tau_membrane(s) işaret değiştirmediği aralıklarda:
  P(s) = sigma_zz(s)^2 + 3 * (tau_surface(s))^2 en fazla 4. dereceden (quartic) bir polinomdur.
- dP/ds = 0 denklemi en fazla 3. dereceden (cubic) türev kökleridir.
- Grid search kesinlikle yasaktır; analitik kök çözümü kullanılır, sanal gürültüler (abs(imag) < 1e-12) filtrelenir.
- Aday noktalar: s = 0, s = L_e, dP/ds = 0 reel kökleri ve tau_membrane(s) = 0 kökleri.
- Global peak stress (max_sigma_vm), konumu (segment_id, s_peak) tam makine hassasiyetinde bulunur.
- Kritik Yük Faktörü: Eğer girdide sigma_yield [F/L^2] sağlanmışsa:
  load_factor = sigma_yield / max_sigma_vm. (sigma_yield verilmemişse None).

## 3. Sayısal Zırhlama ve IEEE-754 Sözleşmesi
- S_omega, B*omega*/Cw ve M_omega*S_omega/Cw işlemlerinde frexp/ldexp mantissa ayrıştırması esastır.
- Taşan değerlerde asla DBL_MAX doyurması yapılmaz; açıkça `OverflowError("Recovered stress overflows IEEE-754 float64 range.")` fırlatılır.
- Sıfırdan farklı alt-normal (< 5e-324) sonuçlarda `FloatingPointError` fırlatılır; simetri sıfırları korunur.
- 10^12 koordinat ötelemesi ve dönme altında skaler gerilmeler 1e-10 bağıl hata ile korunur.

## 4. Bağımsız Benchmark ve Oracle Şartı
- `tests/test_stress_benchmarks.py` içinde standart binary64 np.linalg.solve kullanılmaz; referanslar fractions.Fraction veya 100 basamak Decimal ile çözülür.
