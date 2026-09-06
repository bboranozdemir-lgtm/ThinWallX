# ThinWallX v1.0 — Production-Quality Technical Release

**Belge türü:** Normatif uygulama ve kabul şartnamesi.  
**Hedef:** v1.0 Nihai Üretim Sürümü; v0.1–v0.8 mekanik ve I/O çekirdeği dondurulmuştur.  
**Durum:** Ayrı faz şartnamesi; tamamlanma veya dondurma onayı değildir.  
**Dil:** Python >= 3.10; açık tip ipuçları zorunludur.

Bu dosya mevcut `ACTIVE_PHASE.md` dosyasını kendiliğinden değiştirmez. Uygulama başlangıcı kullanıcı tarafından ayrıca yetkilendirilir. “Zorunludur”, “reddedilir” ve “yasaktır” ifadeleri normatiftir.

---

## 1. Giriş ve Kapsam

### 1.1 Amaç ve Teslimatlar

ThinWallX projesinin v0.1'den v0.8'e kadar uzanan tüm geliştirme aşamalarında ince cidarlı kesit mekaniğinin analitik omurgası ve temel I/O protokolleri inşa edilmiş ve 520 adet bağımsız test ile dondurulmuştur.

v1.0 fazı, kütüphaneyi son kullanıcı mühendisler, araştırmacılar ve CI/CD hatları için tek başına çalışabilir, dağıtılabilir, tam dokümante edilmiş endüstri standardında bir mühendislik aracına (**Production-Quality Technical Release**) dönüştürür.

Bu fazın 4 ana teslimatı bulunmaktadır:
1. **Uçtan Uca CLI (Komut Satırı Arabirimi - `src/thinwallx/cli.py` & `src/thinwallx/__main__.py`):**
   Harici hiçbir üçüncü taraf CLI kütüphanesine (örn. `click`, `typer`) bağımlı olmadan, standart Python `argparse` altyapısıyla çalışan; `inspect`, `analyze`, `plot` ve `convert-dxf` alt komutlarını sunan zırhlı CLI.
2. **Otomatik Mühendislik Hesap Raporu Sentezi (`src/thinwallx/reporting.py`):**
   GitHub-Flavored Markdown formatında, yapısal mühendislik standartlarına uygun hesap çıktısı (Calculation Sheet) üreten raporlayıcı motor.
3. **Kapsamlı Teknik ve Teorik Dokümantasyon (`docs/`):**
   Matematiksel formülasyonları, işaret sözleşmelerini, varsayımları, kısıtları ve benchmark referanslarını içeren eksiksiz teknik doküman paketi.
4. **Paketleme, Dağıtım ve Sürüm Dondurması:**
   `pyproject.toml` sürümünün `1.0.0` yapılması, `[project.scripts]` konsol komutunun (`thinwallx`) tanımlanması ve standart wheel/sdist paketlenebilirliğinin doğrulanması.

---

### 1.2 Dondurulmuş Davranış (Frozen Baseline)

v1.0 fazında önceki fazların mekanik ve I/O hesaplama modülleri kesinlikle dondurulmuştur:

| Faz | Korunacak Sözleşme ve Fonksiyonlar |
|---|---|
| **v0.1** | `Section`, alan, ağırlık merkezi $C$, $I_x, I_y, I_{xy}$, principal properties |
| **v0.2** | Açık kesit enine kayma akışı $q_V$, statik moment integralleri |
| **v0.3** | Açık kesit kayma merkezi $S=(x_s, y_s)$, tork işaret sözleşmesi |
| **v0.4** | Açık kesit Saint-Venant $J$, normalize sektöriyel koordinat $\omega^*$, $C_w$ |
| **v0.5** | `ClosedSection`, çok hücreli Bredt-Batho $H$ matrisi, uyumluluk akıları |
| **v0.6** | `MixedSection`, Tarjan köprü tespiti, çevrim rankı $n_c$, $J_{total}$, sürekli $\omega^*$ |
| **v0.7** | `AppliedLoads`, `calculate_stresses`, $\sigma_{zz}$, $\tau$, $\sigma_{vm}$ analitik tepe taraması |
| **v0.8** | `serialization.py` (F64 hex, UnitSystem), `dxf.py` (ASCII parser), `plotting.py` (Headless plots) |

- v0.1–v0.8 için yazılmış mevcut **520 test** kesinlikle silinemez, atlanamaz veya toleransları gevşetilemez.
- `python -m pytest -W error` koşusunda sıfır hata ve sıfır uyarı korunmalıdır.
- CLI ve Raporlama modülleri mevcut hesaplama çekirdeğini sadece çağırır; yeni mekanik formül veya alternatif çözücü üretemez.

---

### 1.3 Bağımlılıklar ve Kapsam Dışı

- **Çalışma Zamanı Bağımlılıkları:** Yalnızca Python standart kütüphanesi (`argparse`, `json`, `math`, `pathlib` vb.), `numpy` ve isteğe bağlı çizim yolunda `matplotlib`.
- **Kesinlikle Yasak Bağımlılıklar:** `click`, `typer`, `ezdxf`, `scipy`, `sympy`, `shapely`, `jinja2`, `weasyprint`, `pdfkit` veya herhangi bir harici CLI/şablon motoru runtime bağımlılığı olarak eklenemez.
- **Kapsam Dışı:** Etkileşimli GUI/pencereler, Web/REST API, 3D render, sonlu elemanlar çözücüsü, plastik akma yüzeyleri (v1.0 sonrası için ayrılmıştır), dinamik optimizasyon.

---

### 1.4 Modül ve Dosya Sınırları

| Dosya Yolu | Sorumluluk |
|---|---|
| `src/thinwallx/cli.py` | Argparse tabanlı komut satırı arabirimi ve alt komut yönlendiricisi |
| `src/thinwallx/__main__.py` | `python -m thinwallx` desteği için yürütülebilir giriş noktası |
| `src/thinwallx/reporting.py` | Markdown hesap raporu (Calculation Sheet) sentez motoru |
| `src/thinwallx/__init__.py` | Paket sürümünün `1.0.0` yapılması ve public semboller |
| `pyproject.toml` | Sürüm güncellemesi ve `[project.scripts] thinwallx = "thinwallx.cli:main"` tanımı |
| `docs/THEORY_AND_CONVENTIONS.md` | v0.1–v0.7 tüm mekanik formüller, işaret sözleşmeleri ve kısıtlar |
| `docs/CLI_REFERENCE.md` | CLI kullanım kılavuzu, parametreler, alt komutlar ve örnekler |
| `docs/VERIFICATION_BENCHMARKS.md` | Tüm fazların analitik ve bağımsız oracle doğrulama dökümü |
| `tests/test_cli.py` | CLI alt komutları, argüman doğrulama, çıkış kodları ve piping testleri |
| `tests/test_reporting.py` | Markdown hesap raporu yapısı, tablo bütünlüğü ve görsel link testleri |
| `tests/test_packaging.py` | Dağıtım, sürüm bütünlüğü ve CLI console_scripts testleri |

---

## 2. Komut Satırı Arabirimi (CLI) Mimarisi ve API Sözleşmesi

### 2.1 Giriş Noktaları ve Kullanım

CLI hem paket konsol scripti olarak (`thinwallx`) hem de modül olarak (`python -m thinwallx`) çalıştırılabilmelidir.

```bash
thinwallx [--version] [--help] <subcommand> [options]
```

### 2.2 Standart UNIX Çıkış Kodları (Exit Codes)

Hata durumlarında CLI sessiz çökmemeli; açıklayıcı bir hata mesajını `stderr`'e yazarak aşağıdaki standart çıkış kodlarıyla sonlanmalıdır:

| Çıkış Kodu | Anlamı / Hata Kategorisi | Tetikleyen Durumlar |
|:---:|---|---|
| **0** | `SUCCESS` | Komut başarıyla tamamlandı. |
| **1** | `CLI_USAGE_ERROR` | Geçersiz alt komut, eksik zorunlu bayrak, tanımsız argüman. |
| **2** | `FILE_IO_ERROR` | Girdi dosyası bulunamadı, okuma izni yok, geçersiz JSON/DXF dosya biçimi. |
| **3** | `GEOMETRY_ERROR` | Kesit geometri hatası (sıfır uzunluk, negatif kalınlık, geçersiz koordinat). |
| **4** | `TOPOLOGY_ERROR` | Kesit topoloji hatası (bağlantısız grafik, X-kesişimi, geçersiz hücre). |
| **5** | `SINGULAR_ERROR` | Tekil atalet, tekil $H$ matrisi, $C_w=0$ iken bimoment uygulanması. |
| **6** | `NUMERICAL_ERROR` | Float64 aşımı (`OverflowError`), alt-normal sıfırlanma (`FloatingPointError`). |

---

### 2.3 Alt Komutlar (Subcommands)

#### 2.3.1 `inspect` Alt Komutu
Kesit dosyasını (JSON veya DXF) okur ve temel geometrik ve atalet özelliklerini terminale raporlar.

```bash
thinwallx inspect <input_file> [--dxf-unit UNIT] [--thickness THICK] [--json] [--output OUT]
```
- `<input_file>`: `.json` veya `.dxf` uzantılı dosya. Uzantıdan otomatik format tespiti yapılır; aksi halde hata verilir.
- `--dxf-unit`: DXF girdi uzunluk birimi (varsayılan: `mm`).
- `--thickness`: DXF katmanlarında kalınlık belirtilmemişse kullanılacak varsayılan et kalınlığı.
- `--json`: Özellikleri terminale düz metin yerine ThinWallX standart JSON formatında basar.
- `--output`: Çıktıyı terminal (stdout) yerine belirtilen dosyaya yazar.
- **Raporlanan Değerler:** Düğüm sayısı, segment sayısı, topoloji türü (`open`, `closed`, `mixed`), $A$, centroid $C=(x_c, y_c)$, $I_x, I_y, I_{xy}$, asal eksen açısı $\theta_p$, asal atalet momentleri $I_1, I_2$, Saint-Venant $J$, kayma merkezi $S=(x_s, y_s)$, Vlasov çarpılma sabiti $C_w$.

#### 2.3.2 `analyze` Alt Komutu
Kesit geometrisini ve dış yükleri alarak tam lif gerilme analizini icra eder.

```bash
thinwallx analyze <input_file> --loads <loads.json> [options]
```
- Alternatif Yük Girişi: `--loads <file.json>` yerine veya ek olarak CLI üzerinden skaler yükler verilebilir:
  `--N FLOAT`, `--Vx FLOAT`, `--Vy FLOAT`, `--Mx FLOAT`, `--My FLOAT`, `--Tsv FLOAT`, `--B FLOAT`, `--M-omega FLOAT`, `--sigma-yield FLOAT`.
- `--plots-dir <DIR>`: Belirtilirse, analiz sonuçlarına ait 2D geometri (`geometry.png`), kayma akışı (`shear_flow.png`) ve kritik Von Mises gerilme konturu (`stress_vm.png`) grafiklerini bu dizine kaydeder.
- `--report <report.md>`: Belirtilirse, analize ait kapsamlı mühendislik hesap raporunu (Calculation Sheet) Markdown formatında yazar.
- `--output <results.json>`: Analiz sonuçlarını (`StressRecoveryResult`) ThinWallX v0.8 şemasına tam uyumlu JSON olarak kaydeder.
- `--stdout`: Sonuç JSON çıktısını stdout'a yönlendirir (komut zincirleme / piping için).

#### 2.3.3 `plot` Alt Komutu
Mekanik gerilme analizi çalıştırmaksızın doğrudan kesit geometrisini çizdirir.

```bash
thinwallx plot <input_file> --output <geom.png|svg> [--show-thickness] [--show-axes] [--show-ids] [--show-shear-center]
```
- `--output`: Çıktı dosya yolu (`.png` veya `.svg`).
- Bayraklar v0.8 `GeometryPlotOptions` ayarlarına birebir eşlenir.

#### 2.3.4 `convert-dxf` Alt Komutu
Bir DXF dosyasını doğrulanmış, kayıpsız ThinWallX kesit JSON dosyasına dönüştürür.

```bash
thinwallx convert-dxf <input.dxf> --output <output.json> [--source-unit mm] [--target-unit mm] [--thickness THICK] [--section-kind auto|open|closed|mixed]
```
- v0.8 `read_dxf` ve `write_json` modüllerini kullanarak atomik ve güvenli dönüşüm sağlar.

---

### 2.4 Borulama (Piping / Stdin / Stdout Desteği)

- `<input_file>` yerine `-` verildiğinde, JSON verisi doğrudan standart girdiden (`sys.stdin`) okunur.
- `--output` verilmediğinde veya `-` seçildiğinde sonuçlar `sys.stdout`'a aktarılır.
- Hata mesajları, loglar ve ilerleme bildirimleri ASLA `stdout`'a yazılamaz; daima `sys.stderr` kanalına iletilir. Böylece Unix boru hatlarında (`cat sec.json | thinwallx inspect - | jq .`) JSON bozulması engellenir.

---

## 3. Otomatik Mühendislik Hesap Raporu Sentezi (`src/thinwallx/reporting.py`)

### 3.1 Şablon Formatı ve Tasarım İlkeleri

Hesap raporu (Calculation Sheet), yapısal mühendislik tasarım ofislerinde denetim makamlarına (belediye, yapı denetim, bağımsız kontrolör) sunulacak nitelikte, GitHub-Flavored Markdown formatında üretilir.

- Raporlama motoru üçüncü taraf şablon kütüphanelerine (örn. Jinja2) bağımlı olmaksızın, temiz ve katı formatlayıcı fonksiyonlarla oluşturulur.
- Rapor dosyası atomik olarak yazılır (`serialization._atomic_write`); yarım veya bozuk rapor dosyası bırakılamaz.
- Tüm sayısal değerler tutarlı birim etiketleriyle ($L, F, F/L^2$ vb.) birlikte sunulur.

### 3.2 Rapor Bölümleri ve Standart İçerik

Üretilen `report.md` dosyası aşağıdaki hiyerarşik bölümleri içermelidir:

1. **Rapor Başlığı ve Üst Bilgi:**
   - Proje adı (`ThinWallX v1.0 Calculation Report`), üretim tarihi/saati (ISO 8601), girdi dosya yolu, seçilen birim sistemi.
2. **Kesit Geometrisi ve Model Varsayımları:**
   - Segment sayısı, düğüm sayısı, toplam orta-hat uzunluğu $\sum L_i$, minimum ve maksimum et kalınlığı ($t_{min}, t_{max}$).
   - İnce cidarlı modelleme kabulleri (orta-hat integrasyonu, köşe bindirmelerinin ihmali, uniform kalınlık).
   - Topoloji sınıfı: Açık (Tree), Kapalı (Tek/Çok hücreli) veya Karma (Mixed).
3. **Enine Kesit Karakteristik Özellikleri (Tablo):**
   - Alan $A$, Ağırlık merkezi $C=(x_c, y_c)$.
   - Atalet momentleri $I_x, I_y$, Çarpım atalet momenti $I_{xy}$.
   - Asal eksen açısı $\theta_p$, Asal atalet momentleri $I_1, I_2$.
   - Saint-Venant burulma sabiti $J$.
   - Kayma merkezi $S=(x_s, y_s)$.
   - Vlasov sektöriyel çarpılma sabiti $C_w$.
4. **Uygulanan Dış Yükler:**
   - $N$ (Eksenel Kuvvet), $V_x, V_y$ (Enine Kesme Kuvvetleri).
   - $M_x, M_y$ (Eğilme Momentleri), $T_{sv}$ (Saint-Venant Burulma Momenti).
   - $B$ (Vlasov Bimomenti), $M_\omega$ (İkincil Çarpılma Burulma Momenti).
   - $\sigma_{yield}$ (Tasarım Akma Dayanımı - varsa).
5. **Kritik Gerilme Analizi ve Kapasite Kontrolü:**
   - Maksimum eşdeğer gerilme: $\sigma_{vm, max}$.
   - Tepe noktasının konumu: Segment ID, lokal $s_{peak}$ koordinatı ve global koordinatlar $(x_{peak}, y_{peak})$.
   - Akma Güvenlik Faktörü (Load Factor): $\lambda = \sigma_{yield} / \sigma_{vm, max}$ (akma dayanımı tanımlıysa).
   - Sonuç bileşkeleri doğrulaması tablosu ($N_{rec}, M_{x,rec}, M_{y,rec}, B_{rec}$ ve denge hata farkları).
6. **Grafiksel Çıktılar (Görsel Ekler):**
   - Eğer grafikler üretilmişse göreli linkler ile rapora gömülür:
     * `![Kesit Geometrisi](./geometry.png)`
     * `![Kayma Akisi Dagilimi](./shear_flow.png)`
     * `![Eşdeger Gerilme Dagilimi](./stress_vm.png)`

---

## 4. Dokümantasyon Şartnamesi (`docs/`)

v1.0 sürümü ile birlikte kod tabanında 3 adet temel başvuru dokümanı tamamlanmış ve dondurulmuş olarak bulunmalıdır:

### 4.1 `docs/THEORY_AND_CONVENTIONS.md`
- İnce cidarlı kesit teorisinin temel diferansiyel ve integral denklemleri.
- 2D kartezyen koordinat sistemi (sağ el kuralı: $x$ sağa, $y$ yukarı, $z$ kesit dışına doğru).
- Pozitif iç kuvvet ve moment tanımları ($N, V_x, V_y, M_x, M_y, T_{sv}, B, M_\omega$).
- Düz segment kapalı formül integrasyon matrisleri (alan, moment, atalet).
- Bredt-Batho hücre dolaşım uyumluluğu ($H \boldsymbol{\phi} = 2 \mathbf{A}$).
- Tarjan köprü ayrıştırması, biconnected bloklar ve çevrim rankı teoremi ($n_c = |E_c| - |V_c| + k_c$).
- Vlasov çarpılma teorisi, sektöriyel alan sürekliliği ve sıfır-ortalama şartı ($\int \omega^* dA = 0$).
- Kuartik tepe gerilme polinomu ve dP/ds kök arama kuralları.
- Sayısal zırhlama stratejisi (`frexp`/`ldexp` mantis ayrıştırması, IEEE-754 sınırları).

### 4.2 `docs/CLI_REFERENCE.md`
- CLI kurulumu, genel komut yapısı ve alt komutların ayrıntılı parametre listesi.
- Çıkış kodları (Exit codes) ve hata yakalama senaryoları.
- Borulama (piping) ve otomasyon senaryoları (örnek Bash/PowerShell betikleri).
- Standart DXF ve JSON girdi formatı örnekleri.

### 4.3 `docs/VERIFICATION_BENCHMARKS.md`
- v0.1–v0.8 boyunca uygulanan tüm bağımsız test oracle'larının özeti:
  * Analitik açık I-kesiti, L-köşebent, kanal ve U-profilleri.
  * Çok hücreli kapalı kutu kirişler (Box girders).
  * Karma köprü-hücre kombinasyonları (Barbell, şapka/omega profilleri).
  * Kıyaslama toleransları ($rtol \le 10^{-10}$–$10^{-12}$) ve literatür referansları.

---

## 5. Paketleme ve Dağıtım Şartnamesi

1. **`pyproject.toml` Yapılandırması:**
   - Sürüm numarası tam olarak `1.0.0` yapılmalıdır:
     ```toml
     [project]
     name = "thinwallx"
     version = "1.0.0"
     ```
   - Konsol betiği giriş noktası eklenmelidir:
     ```toml
     [project.scripts]
     thinwallx = "thinwallx.cli:main"
     ```
   - İsteğe bağlı bağımlılıklar (`plots`, `test`) korunmalıdır.
2. **Paketlenebilirlik Güvencesi:**
   - Standart Python paketleme araçlarıyla (`python -m build` veya `pip install -e .`) hatasız wheel ve sdist arşivleri üretilebilmelidir.
   - Paketin kurulmasının ardından `thinwallx --version` çağrısı `ThinWallX 1.0.0` metnini döndürmelidir.

---

## 6. Doğrulama ve Test Matrisi (T01 – T40)

Yeni testler `tests/test_cli.py`, `tests/test_reporting.py` ve `tests/test_packaging.py` dosyalarında yer almalıdır.

| Test ID | Kategori | Doğrulanacak Kural / Senaryo |
|:---:|---|---|
| **T01** | CLI Core | `thinwallx --version` ve `thinwallx --help` doğru metin ve exit code 0 döner. |
| **T02** | CLI Core | Tanımsız alt komut veya geçersiz bayrak `exit code 1` üretir. |
| **T03** | CLI Core | Olmayan dosya yolu verildiğinde `exit code 2` üretir. |
| **T04** | CLI Inspect | JSON açık kesit okuma: $A, I_x, I_y, S, J, C_w$ terminale doğru basılır. |
| **T05** | CLI Inspect | DXF kapalı kutu okuma: Çok hücreli özellikler doğru hesaplanır. |
| **T06** | CLI Inspect | `--json` bayrağı ile çıktının standart JSON olması doğrulanır. |
| **T07** | CLI Inspect | Geçersiz/bozuk geometrili kesitte `exit code 3` (GeometryError) üretir. |
| **T08** | CLI Analyze | JSON kesit + JSON yükler: Gerilme analizi hatasız tamamlanır (`exit code 0`). |
| **T09** | CLI Analyze | DXF kesit + CLI yük argümanları (`--N`, `--Mx`, vb.) başarıyla çözülür. |
| **T10** | CLI Analyze | `--plots-dir` ile 3 grafiğin (geom, shear, stress) dizine kaydedildiği doğrulanır. |
| **T11** | CLI Analyze | `--report` ile üretilen Markdown raporunun diskte oluştuğu doğrulanır. |
| **T12** | CLI Analyze | `--output` ile sonuçların geçerli ThinWallX JSON formatında yazıldığı doğrulanır. |
| **T13** | CLI Analyze | $C_w=0$ kesite bimoment $B$ verildiğinde `exit code 5` (SingularError) üretir. |
| **T14** | CLI Analyze | Float64 aşımı yaratan aşırı yükte `exit code 6` (NumericalError) üretir. |
| **T15** | CLI Convert | DXF dosyasını hatasız JSON kesitine dönüştürür ve geri okunabilirliğini sınar. |
| **T16** | CLI Plot | Sadece geometri çizimi komutu başarıyla PNG ve SVG üretir. |
| **T17** | CLI Piping | `stdin` üzerinden JSON kesit besleme (`thinwallx inspect -`) doğrulanır. |
| **T18** | CLI Piping | `--stdout` ile analiz sonuçlarının stdout'a temiz aktarıldığı doğrulanır. |
| **T19** | CLI Safety | Matplotlib kurulu değilken `inspect` ve `analyze` (plotsuz) çalışır; `plot` açıklayıcı hata verir. |
| **T20** | CLI Safety | Çıktı dosyasında `overwrite=False` iken üzerine yazma engellenir. |
| **T21** | Reporting | Rapor başlık, tarih, girdi yolu ve birim meta verilerini eksiksiz içerir. |
| **T22** | Reporting | Geometri ve modelleme kabulleri bölümü segment sayılarını doğru yansıtır. |
| **T23** | Reporting | Enine kesit özellikleri tablosu $A, C, I, J, S, C_w$ değerlerini tam hassasiyetle yansıtır. |
| **T24** | Reporting | Uygulanan yükler tablosu tüm yük bileşenlerini doğru birimlerle listeler. |
| **T25** | Reporting | Kritik gerilmeler tablosu tepe $\sigma_{vm}$, konum ve akma faktörünü tam yansıtır. |
| **T26** | Reporting | Sonuç bileşkeleri denge kontrol tablosu ($N, M_x, M_y, B$) doğrulanır. |
| **T27** | Reporting | Üretilen grafiklerin Markdown resim linkleri (`![...](./...)`) geçerli ve göreli yolludur. |
| **T28** | Reporting | Rapor oluşturma atomik yazma kullanır; işlem başarısızlığında yarım dosya kalmaz. |
| **T29** | Reporting | Asimetrik kesit ($I_{xy} \ne 0$) için rapor tablosunda asal açılar ve ataletler doğrudur. |
| **T30** | Reporting | Çok hücreli kapalı ve karma kesitler için rapor eksiksiz üretilir. |
| **T31** | Packaging | `pyproject.toml` sürümü `1.0.0` ile `thinwallx.__version__` birebir eşleşir. |
| **T32** | Packaging | `thinwallx` console script giriş noktası `thinwallx.cli:main` olarak tanımlıdır. |
| **T33** | Packaging | `python -m thinwallx` doğrudan çalıştırılabilir. |
| **T34** | Packaging | Kaynak kodlarda çözülmemiş `TODO`, `FIXME` veya placeholder bulunmaz. |
| **T35** | Packaging | `docs/` altındaki 3 ana dokümanın tam, eksiksiz ve bağlantılarının sağlam olduğu doğrulanır. |
| **T36** | Packaging | Dağıtım arşivi oluşturma simülasyonu başarıyla geçer. |
| **T37** | Regression | v0.1 Temel geometri ve atalet testleri (tamamı geçer). |
| **T38** | Regression | v0.2–v0.4 Açık kesit kayma akışı, kayma merkezi ve çarpılma testleri (tamamı geçer). |
| **T39** | Regression | v0.5–v0.6 Kapalı ve karma kesit mekaniği testleri (tamamı geçer). |
| **T40** | Regression | v0.7–v0.8 Bileşik gerilme, JSON, DXF ve Çizim testleri (tamamı geçer, toplam 520+ test). |

---

## 7. Kabul Kriterleri (Acceptance Criteria Checklist)

v1.0 sürümünün resmî olarak dondurulabilmesi için aşağıdaki tüm maddelerin istisnasız sağlanması zorunludur:

- [ ] `pyproject.toml` sürümü `1.0.0` yapılmış ve `console_scripts` tanımlanmıştır.
- [ ] `src/thinwallx/cli.py` bağımsız `argparse` ile `inspect`, `analyze`, `plot`, `convert-dxf` komutlarını hatasız sunmaktadır.
- [ ] Standart UNIX çıkış kodları (0–6) tam olarak uygulanmış ve test edilmiştir.
- [ ] `src/thinwallx/reporting.py` profesyonel GitHub-Flavored Markdown hesap raporunu eksiksiz üretmektedir.
- [ ] `docs/THEORY_AND_CONVENTIONS.md`, `docs/CLI_REFERENCE.md` ve `docs/VERIFICATION_BENCHMARKS.md` eksiksiz oluşturulmuştur.
- [ ] Önceki fazlara ait dondurulmuş kodlar (v0.1–v0.8) değiştirilmemiş ve geriye dönük tam uyumluluk korunmuştur.
- [ ] Sıfır ek harici kütüphane kuralı korunmuştur (yalnızca `numpy`, `pytest` ve isteğe bağlı `matplotlib`).
- [ ] `tests/test_cli.py`, `tests/test_reporting.py` ve `tests/test_packaging.py` test paketleri yazılmıştır.
- [ ] `python -m pytest -W error` komutu, tüm eski 520 test ve yeni v1.0 testleri dahil olmak üzere **sıfır hata ve sıfır uyarı** ile başarıyla tamamlanmaktadır.
- [ ] Faz tamamlama raporu sunulmuş ve **STOP** kuralına uyulmuştur.
