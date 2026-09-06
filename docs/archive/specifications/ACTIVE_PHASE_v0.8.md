# ThinWallX v0.8 — JSON, ASCII DXF Orta-Hat İçe Aktarımı ve Teknik Görselleştirme

**Belge türü:** Normatif uygulama ve kabul şartnamesi.  
**Hedef:** v0.8; v0.1–v0.7 mekanik çekirdeği dondurulmuştur.  
**Durum:** FROZEN — kullanıcı tarafından resmen onaylandı; paket sürümü 0.8.0.  
**Dil:** Python >= 3.10; açık tip ipuçları zorunludur.

Bu dosya mevcut ACTIVE_PHASE.md dosyasını değiştirmez. Uygulama başlangıcı kullanıcı tarafından ayrıca yetkilendirilir. “Zorunludur”, “reddedilir” ve “yasaktır” ifadeleri normatiftir.

## 1. Giriş ve Kapsam

### 1.1 Amaç ve teslimatlar

v0.8, doğrulanmış ince cidarlı kesit mekaniğini dış veri kaynaklarına ve statik mühendislik çizimlerine bağlar:

1. Section, ClosedSection, MixedSection, AppliedLoads ve StressRecoveryResult için sürümlü, katı ve kayıpsız JSON.
2. ASCII DXF orta-hatlarından doğrulanmış açık, kapalı veya karma kesit elde edilmesi.
3. GUI açmadan geometri, kayma akışı ve gerilme dağılımlarının PNG/SVG çıktıları.

Bu üç modül adaptördür; yeni mekanik çözücü değildir. A, C, Ix, Iy, Ixy, asal eksenler, q, S, J, omega*, Cw ve gerilmeler mevcut doğrulanmış API'lerden alınır.

### 1.2 Dondurulmuş davranış

| Faz | Korunacak sözleşme |
|---|---|
| v0.1 | Düz segment integralleri, alan, centroid, atalet ve atan2 asal eksenleri |
| v0.2 | Açık kesit enine kayma akışı |
| v0.3 | Açık kesit kayma merkezi ve tork işaretleri |
| v0.4 | Saint-Venant J, normalize sektöriyel koordinat, Cw |
| v0.5 | Kapalı çok hücreli Bredt–Batho ve uyumluluk |
| v0.6 | Karma topoloji, bağımsız hücre blokları, J_total, sürekli omega* |
| v0.7 | Yükler, bileşik gerilmeler, analitik maksimum ve sonuç bileşkeleri |

- Hesap modülleri, mevcut toleranslar, istisnalar, yapıcı varsayılanları ve test beklentileri değiştirilmez.
- Yeni API'ler ayrı modüllerde yer alır. Monkey-patch veya eski sınıflara dinamik metot enjeksiyonu yoktur.
- Paket dışa aktarımları ve isteğe bağlı grafik bağımlılığı bildirimi gerekirse yalnız eklemeli değişiklik olabilir.
- Kullanıcının bildirdiği taban 389 başarılı testtir. Uygulama başlangıcında gerçek taban yeniden ölçülür; test sayısını sağlamak için test silmek veya atlamak yasaktır.
- Çekirdekte yeni hata bulunursa adaptörde gizlenmez; ayrı raporlanır ve yetkisiz mekanik düzeltme yapılmaz.
- Bu şartname dosyası v0.8 uygulamasının başlatıldığı anlamına gelmez.

### 1.3 Bağımlılıklar ve kapsam dışı

Çalışma zamanı Python standart kütüphanesi, numpy ve yalnız çizim yolunda matplotlib ile sınırlıdır. Testler pytest, fractions ve decimal kullanabilir.

ezdxf, sympy, scipy, shapely, CAD motorları, FEM/ağ kütüphaneleri ve üçüncü taraf JSON doğrulama çalışma zamanı eklenemez. Matplotlib'in kendi dağıtım bağımlılıkları dışında yeni doğrudan grafik bağımlılığı yoktur.

DXF dışa aktarımı, binary DXF, DWG, INSERT/blok açılımı, yay yaklaşımı, spline örneklemesi, 3D projeksiyon, değişken et kalınlığı, katı model, plastik analiz, burkulma, optimizasyon, GUI ve web arayüzü kapsam dışıdır.

Çizim için analitik alan örneklenebilir; örneklerden rijitlik, gerilme katsayısı veya mekanik maksimum hesaplanamaz.

### 1.4 Modül ve teslimat sınırları

| Yol | Sorumluluk |
|---|---|
| src/thinwallx/serialization.py | Katı sözlük/JSON codec'i |
| src/thinwallx/dxf.py | Grup çifti okuyucu, entity ayrıştırıcı, normalizasyon |
| src/thinwallx/plotting.py | Tembel matplotlib yükleme ve statik grafikler |
| schemas/thinwallx-0.8.schema.json | Bu belgenin yapısal JSON Schema karşılığı |
| tests/test_serialization.py | Round-trip ve güvenilmeyen girdi testleri |
| tests/test_dxf_import.py | Sözdizimi, geometri ve topoloji |
| tests/test_dxf_benchmarks.py | Bağımsız mekanik oracle'lar |
| tests/test_plotting.py | Headless çıktı, determinizm, kaynak sahipliği |
| tests/fixtures/dxf/ | Küçük, incelenebilir ASCII fixture'lar |

Yeni yardımcı tipler küçük dataclass'lardır. Gelecek fazlar için eklenti mimarisi veya genel sınıf fabrikası kurulmaz.

### 1.5 Birimler

Çekirdek birimsiz değil, birimden bağımsızdır. Girdiler tutarlı L uzunluk ve F kuvvet birimleri taşır.

| Büyüklük | Boyut |
|---|---|
| x, y, s, t, node_tolerance | L |
| A, omega* | L² |
| Ix, Iy, Ixy, J | L⁴ |
| Cw | L⁶ |
| N, Vx, Vy | F |
| Mx, My, Tsv, M_omega | F L |
| B | F L² |
| q | F/L |
| sigma, tau, sigma_yield | F/L² |

Sağ elli x–y–z sistemi ve p1→p2 segment yönü korunur. Birim dönüşümü ile koordinat dönüşümü, yeniden etiketleme işlemiyle karıştırılmaz.

## 2. JSON Şeması ve API Sözleşmesi

### 2.1 Format

UTF-8, BOM'suz, tek JSON kök nesnesi kullanılır. Kök alanların tamamı zorunludur:

| Alan | Tür/değer |
|---|---|
| format | Tam olarak "thinwallx" |
| schema_version | Tam olarak "0.8" |
| number_encoding | Tam olarak "float64-hex" |
| kind | "section", "applied_loads" veya "stress_result" |
| units | §2.4 birim kaydı |
| data | kind'a bağlı kapalı kayıt |

Tüm kayıtlar kapalı şemadır: additionalProperties=false. Bilinmeyen alan, eksik zorunlu alan, yinelenen anahtar, yanlış tür, bilinmeyen sürüm ve enum reddedilir. Gelecek sürümler tahmin edilmez.

Yapısal JSON Schema dosyası aynı kuralları açıklar. Çapraz referans, sonluluk ve topoloji ayrıca Python'da denetlenir; üçüncü taraf validator zorunlu değildir.

Standart json okuyucusunda yinelenen anahtarlar object_pairs_hook, standart dışı sayılar parse_constant ile reddedilir. Yazıcı allow_nan=False kullanır. Python'un JSON anahtarlarını string'e çevirme davranışı ID serileştirmesinde kullanılmaz. [Python JSON belgeleri](https://docs.python.org/3/library/json.html)

### 2.2 F64: bit düzeyinde kayıpsız fiziksel sayılar

Normatif fiziksel sayı kaydı:

~~~json
{"f64":"0x1.0000000000000p+0"}
~~~

Metin sonlu binary64 değerin float.hex() çıktısıdır. Okuyucu float.fromhex() sonrasında tekrar hex() alarak kanonik metin eşitliğini doğrular. Böylece gereksiz basamaklı veya aralık dışı olup sessiz yuvarlanan hex girdisi reddedilir.

| Değer | Kanonik metin |
|---|---|
| +0 | 0x0.0p+0 |
| -0 | -0x0.0p+0 |
| 1 | 0x1.0000000000000p+0 |
| En küçük pozitif alt-normal | 0x0.0000000000001p-1022 |
| En büyük sonlu binary64 | 0x1.fffffffffffffp+1023 |

Fiziksel alanlarda çıplak JSON number veya serbest decimal string kabul edilmez. Tek normatif kodlama vardır. Amaç JSON tüketicisinin sayı türünden bağımsız hassasiyettir. Hex dönüşümü binary64 değerini tam taşır. [Python float.hex/fromhex](https://docs.python.org/3/library/stdtypes.html#float.hex)

\[
\operatorname{pack}_{64}(x_{\rm out})=\operatorname{pack}_{64}(x_{\rm in}).
\]

Normal sayılarda bu koşul istenen \(10^{-15}\) hata sınırından güçlüdür. Sıfır ve alt-normalde bağıl hata yerine bit eşitliği uygulanır. Signed zero korunur. Codec içinde fiziksel ölçekleme, yuvarlama veya “küçük değer temizleme” yoktur.

NaN ve her türlü sayı-sonsuzluk literal'i reddedilir. Tek semantik istisna, sıfır gerilmede v0.7 yük faktörüdür:

~~~json
{"special":"positive_infinity"}
~~~

Bu kayıt yalnız stress_result.data.load_factor alanında kullanılabilir; max_sigma_vm sıfır olmalıdır. Diğer alanlara taşan değerleri bu etiketle saklamak yasaktır. Bu bir JSON sayı sabiti değildir. [RFC 8259](https://www.rfc-editor.org/rfc/rfc8259)

### 2.3 PortableId

Mevcut ID alanları Any olsa da her Python nesnesinin taşınabilir JSON karşılığı olduğu iddia edilmez. Desteklenen küme:

| Python tipi | Wire kayıt |
|---|---|
| None | {"type":"none"} |
| str | {"type":"str","value":"web-A"} |
| int; bool hariç | {"type":"int","value":"17"} |
| bool | {"type":"bool","value":true} |
| Sonlu float | {"type":"float","value":{"f64":"0x1.0000000000000p+0"}} |
| tuple | {"type":"tuple","items":[{"type":"str","value":"A"},{"type":"int","value":"2"}]} |

Integer metni kanonik decimal'dir: "0" veya isteğe bağlı eksi ve sıfırla başlamayan basamaklar; +1, -0 ve 01 kabul edilmez. Tuple özyinelemeli aynı kümeyi kullanır.

Özel sınıf, bytes, set, list/dict ID ve sonlu olmayan float TypeError üretir. repr ile kimlik türetme, pickle, eval, dinamik import yoktur.

Node/Segment ID'si topolojik referans değildir. Ref alanı ayrıdır. Aynı ID'li ayrı uçlar korunur. Kayıpsızlık desteklenen tür, kamusal değer, sıra ve yön içindir; bellek adresi, aliasing ve tembel önbellekler için değildir.

Result sözlüğünde Python eşitliği bakımından çakışan anahtarlar (1, True, 1.0) reddedilir. Wire etiketlerinin farklı olması tek başına yeterli değildir.

### 2.4 Birim kaydı

units yalnız şu iki zorunlu alanı taşır:

- length: "m", "mm", "cm", "in", "ft" veya "unspecified".
- force: "N", "kN", "lbf" veya "unspecified".

Etiketler metadata'dır; codec sayı dönüştürmez. unspecified, SI varsayımı değildir. Geçersiz yazım reddedilir.

Çekirdek nesneleri birim saklamadığından çözücü bir DecodedDocument döndürür: value ve units. Bu zarf sayesinde dosya yeniden yazılırken birimler kaybolmaz.

### 2.5 Kesit kaydı

kind="section" data alanları:

| Alan | Kural |
|---|---|
| class | "Section", "ClosedSection" veya "MixedSection" |
| node_tolerance | F64, sonlu ve >=0, mutlak uzunluk |
| safety_factor | Section için null; diğer sınıflarda mevcut yapıcıya uygun F64 |
| nodes | En az iki NodeRecord |
| segments | En az bir SegmentRecord |

NodeRecord tam alanları: ref:string, id:PortableId, x:F64, y:F64.  
SegmentRecord tam alanları: id:PortableId, p1:string, p2:string, t:F64.

- Ref'ler benzersiz ve boş olmayan string'dir. Tüm p1/p2 tanımlı olmalıdır.
- Kullanılmayan düğüm kaydı, kendine bağlı kenar, sıfır/sonsuz uzunluk ve pozitif olmayan t reddedilir.
- Segment sırası ve uç yönü korunur.
- Yazıcı her segmentin iki ucunu özgün koordinatlarıyla sıra içinde n0,n1,n2 şeklinde kaydeder. Ortak nesneleri ayrı kayıtlarla yazması serbest değil, kanonik yazıcı davranışıdır. Geometrik küme temsilcisi ham koordinatın yerine konmaz.
- Okuyucu paylaşılan ref'leri kabul eder; kanonik yeniden yazmada ayrı endpoint kayıtlarına açılabilir. Metin değil, nesne değerleri round-trip kapsamındadır.
- Okuyucu mevcut sınıfı validate=True ile kurar. Wire kaydında validate alanı yoktur.
- Yazıcı da kaynak kesidi doğrular; validate=False ile yaratılmış geçersiz nesne arşivlenmez.
- MixedSection saf açık grafik taşısa bile aynı sınıfta korunur. Sınıf JSON okurken otomatik seçilmez.
- H, H_scaled, J, Cw, centroid ve önbellekler kaydedilmez ve yazmak için tetiklenmez. Bunlar geometriden türetilir.
- JSON okuyucusu snapping, T-splitting veya onarım yapmaz.

### 2.6 AppliedLoads kaydı

Tüm yükler açıkça zorunludur; eksik alan sıfıra tamamlanmaz.

| Alan | Tip | Boyut |
|---|---|---|
| N, Vx, Vy | Her biri F64 | F |
| Mx, My, Tsv, M_omega | Her biri F64 | F L |
| B | F64 | F L² |
| sigma_yield | null veya pozitif F64 | F/L² |

AppliedLoads doğrulaması uygulanır; momentler yeniden merkezlenmez.

Tam bir geçerli belge:

~~~json
{
  "format":"thinwallx",
  "schema_version":"0.8",
  "number_encoding":"float64-hex",
  "kind":"applied_loads",
  "units":{"length":"mm","force":"N"},
  "data":{
    "N":{"f64":"0x0.0p+0"},
    "Vx":{"f64":"0x0.0p+0"},
    "Vy":{"f64":"0x0.0p+0"},
    "Mx":{"f64":"0x0.0p+0"},
    "My":{"f64":"0x0.0p+0"},
    "Tsv":{"f64":"0x0.0p+0"},
    "B":{"f64":"0x0.0p+0"},
    "M_omega":{"f64":"0x0.0p+0"},
    "sigma_yield":null
  }
}
~~~

### 2.7 StressRecoveryResult kaydı

kind="stress_result" data alanları eksiksiz olarak:

| Alan | Tip/kural |
|---|---|
| segment_profiles | Boş olmayan sıralı ProfileEntry dizisi |
| max_sigma_vm | F64, >=0 |
| peak_segment_id | Mevcut profile anahtarı olan PortableId |
| peak_s | F64, seçili segmentin [0,L] aralığı |
| peak_location_xy | Tam iki F64 |
| load_factor | null, pozitif F64 veya §2.2 sonsuzluk etiketi |
| resultant_N | F64 |
| resultant_Mx | F64 |
| resultant_My | F64 |
| resultant_B | F64 |

ProfileEntry tam alanları: key:PortableId, profile:ProfileRecord.

ProfileRecord tam alanları:

| Alan | Tip/anlam |
|---|---|
| segment_id | PortableId; key ile tür/değer eşitliği |
| length | Pozitif F64 |
| thickness | Pozitif F64 |
| sigma_zz_coeffs | İki F64; [a,b], sigma(s)=a s+b |
| tau_membrane_coeffs | Üç F64; [a,b,c], tau(s)=a s²+b s+c |
| tau_sv_surface | F64, >=0 |
| peak_sigma_vm | F64, >=0 |
| s_peak | F64, [0,length] |

Katsayılar xi değil boyutlu s içindir. Normal katsayıların boyutları F/L³ ve F/L²; membran katsayılarının F/L⁴, F/L³ ve F/L²'dir. s, özgün p1→p2 yönünde ölçülür.

ProfileEntry sırası mevcut sonuç sözlüğünün insertion sırasını korur. Anahtarlar JSON object key'ine çevrilmez.

Çözücü kuralları:

1. SegmentStressProfile mevcut kurucuyla oluşturulur. Özel _sigma_xi ve _tau_xi dosyadan yüklenmez.
2. Yeniden hesaplanan peak_sigma_vm ve s_peak, arşivle aynı desteklenen hesap ortamında bit eşit olmalıdır. Değilse ValueError; frozen cache alanına zorla yazma yoktur.
3. Farklı nümerik sürümde yeniden hesap farklıysa açık uyumluluk hatası verilir. Her platform/sürüm arasında bit eşitliği vaat edilmez.
4. Global maksimum profile maksimumlarından denetlenir. Eşit maksimumlar varsa kaydedilmiş seçili profile korunabilir; o profile global maksimumu taşımalıdır.
5. peak_s seçili profile s_peak ile eşleşmelidir.
6. Sonuç nesnesinin içinde geometri veya yük bulunmaz. Standalone kayıt peak_location_xy ve sonuç bileşkelerinin mekanik doğruluğunu ispatlamaz; bunları kayıpsız saklar.
7. Kesitle çizimde ID, uzunluk, t ve tepe koordinatı eşleşmesi ayrıca denetlenir. Yük faktörü–akma gerilmesi ilişkisi ancak yükler ayrıca varsa doğrulanabilir.
8. Yapısal dosya doğrulaması bağımsız mekanik denetim veya kriptografik özgünlük değildir.

### 2.8 API ve dosya sahipliği

Aşağıdakiler imza sözleşmesidir; bu belge uygulama gövdesi içermez.

~~~text
Serializable = Section | ClosedSection | MixedSection | AppliedLoads | StressRecoveryResult
to_dict(value: Serializable | DecodedDocument, *, units: UnitSystem | None = None) -> dict[str, object]
from_dict(data: Mapping[str, object], *, limits: DecodeLimits | None = None) -> DecodedDocument
to_json(value: Serializable | DecodedDocument, *, units: UnitSystem | None = None, indent: int | None = None) -> str
from_json(text: str, *, limits: DecodeLimits | None = None) -> DecodedDocument
write_json(value: Serializable | DecodedDocument, path: str | Path, *, units: UnitSystem | None = None, overwrite: bool = False) -> Path
read_json(path: str | Path, *, limits: DecodeLimits | None = None) -> DecodedDocument
~~~

UnitSystem frozen dataclass alanları length ve force; ikisi de varsayılan unspecified. DecodedDocument frozen dataclass alanları value ve units'tir.

Ham nesnede units verilmezse unspecified kullanılır. Zarf yeniden yazılırken kendi birimi korunur; farklı units argümanı sessiz yeniden etiketleme olacağından reddedilir.

Modül fonksiyonları istenen to/from sözleşmesini sağlar; eski sınıflara metot eklemek gerekmez.

- to_dict bağımsız kayıt üretir; kayıt mutasyonu kaynağı değiştirmez. Gereksiz tam deepcopy yapılmaz.
- from_dict girdiyi değiştirmez. Bool fiziksel sayı gibi kabul edilmez.
- Canonical JSON: sort_keys=True, compact separators, ensure_ascii=False, indent=None. Diziler sıralanmaz.
- Unicode etiketleri normalize edilmez; geçersiz tekil surrogate reddedilir.
- to_json dosya yazmaz. write_json UTF-8, tek son LF ve varsayılan overwrite=False kullanır.
- Dosya yazımı aynı dizindeki geçici dosyadan atomik tamamlanır. Hata halinde hedefte kısmi belge bırakılmaz.
- İçerikten dosya yolu, URL veya yürütülecek komut çözülmez.

### 2.9 Hatalar ve sınırlar

| Durum | Hata |
|---|---|
| Desteklenmeyen Python tipi | TypeError |
| JSON sözdizimi/şema/ref/sürüm | ValueError; JSONDecodeError korunabilir |
| Sayısal geçersiz geometri | GeometryError |
| Geçersiz topoloji | Mevcut TopologyError davranışı |
| Temsil aralığı | OverflowError veya FloatingPointError |
| Dosya erişimi | Özgün OSError alt türü |

DecodeLimits pozitif tam sayılı alanları: max_bytes=16777216, max_depth=64, max_nodes=10000, max_segments=5000, max_string_chars=4096, max_profiles=5000.

Byte sınırı parse öncesi; diğer sınırlar nesne kurmadan önce uygulanır. Derinlik sınırı için string içindeki parantezleri saymayan, kaçışları bilen ön tarama yapılmalıdır; yalnız parse sonrası kontrol bellek/stack koruması sayılmaz. Sınır aşımı açıklayıcı ValueError'dır. Kullanıcı sonlu limitleri açıkça artırabilir.

Hatalar data.segments[3].t gibi alan yolunu içerir; bütün dosya loglanmaz.

## 3. ASCII DXF Parser Matematiksel ve Geometrik Temelleri

### 3.1 Desteklenen format

$ACADVER zorunludur; AC1009 ve AC1015 desteklenir. Binary DXF imzası reddedilir.

| Entity | AC1009 | AC1015 |
|---|---|---|
| LINE | Evet | Evet |
| 2D POLYLINE/VERTEX/SEQEND | Evet | Evet |
| LWPOLYLINE | Hayır | Evet |

Seçili model-space katmanındaki ARC, CIRCLE, ELLIPSE, SPLINE, INSERT, HATCH, REGION, SOLID, 3DFACE ve diğer desteklenmeyen entity'ler reddedilir. Yaylar kirişe çevrilmez.

HEADER okunur. TABLES, BLOCKS, CLASSES ve OBJECTS yalnız yapısal geçilir; blok genişletilmez. ENTITIES içindeki INSERT sessiz atlanamaz.

### 3.2 Durum makinesi ve güvenli parsing

- Dosya grup-kod/değer satır çiftleriyle okunur. Tek kalan satır, geçersiz kod ve eksik değer hatadır.
- Grup kodu decimal integer; sayısal değer ASCII noktalı decimal ve isteğe bağlı bilimsel üstür. Virgüllü decimal, NaN ve inf reddedilir.
- SECTION/section-name/ENDSEC/EOF sıralaması doğrulanır. İç içe section, ikinci ENTITIES, yinelenen HEADER, eksik EOF ve EOF sonrası veri reddedilir; son boş satırlar kabul edilir.
- POLYLINE alt durumu VERTEX dizisini SEQEND ile kapatır. Orphan VERTEX/SEQEND ve yarım zincir reddedilir.
- Tekil geometrik grup kodu tekrarı hatadır. LWPOLYLINE tekrar eden 10/20 kayıtları ayrı vertex'tir; son değerle ezilmez.
- Grup 90 sayısı gerçek vertex sayısıyla eşleşir. Açık zincir en az iki, kapalı en az üç ayrı köşe taşır.
- 999 açıklamaları geometriye etki etmez.
- Entity metadata allowlist'i: 5,6,8,48,60,62,67,100,102 dengeli gruplar,330,370,410,420,430,440,999. İlgili entity'nin aşağıdaki geometrik kodları ayrıca desteklenir. Bilinmeyen geometrik kod, XDATA ve extension dictionary bu dar profilde reddedilir.
- Handle varsa geçerli hexadecimal ve dosyada tekildir. Entity eşlemesi için kullanılır; canonical fiziksel ID değildir.
- Hata satır numarası, entity türü, handle ve layer bağlamını içerir.

Bu dar profil bütün CAD dosyalarını kabul etmeyi vaat etmez. Açık reddetme, sessiz geometri kaybına tercih edilir.

### 3.3 DXF grup kodları

| Kod | LINE | LWPOLYLINE | POLYLINE / VERTEX |
|---|---|---|---|
| 0 | Tür | Tür | Başlık/vertex/son |
| 5 | Handle | Handle | Handle |
| 8 | Layer | Layer | Ana layer |
| 10,20 | p1 x,y | Sıralı vertex x,y | Vertex x,y; başlık dummy x,y=0 |
| 30 | p1 z; yoksa 0 | Desteklenmez | Vertex z / başlık elevation |
| 11,21,31 | p2 x,y,z; z yoksa 0 | Geçersiz | Geçersiz |
| 38 | Geçersiz | Ortak elevation; yoksa 0 | Geçersiz |
| 39 | Extrusion thickness; yalnız 0 | Yalnız 0 | Yalnız 0 |
| 40,41,43 | Geçersiz | Width; yalnız 0 | İlgili default/vertex width; yalnız 0 |
| 42 | Geçersiz | Bulge; yalnız 0 | Vertex bulge; yalnız 0 |
| 66 | Geçersiz | Geçersiz | Eski entities-follow; SEQEND yerine geçmez |
| 70 | Geçersiz | Bit 1 closed,128 plinegen | Başlık bit 1,128; diğer geometrik bitler reddedilir |
| 90 | Geçersiz | Vertex sayısı | Geçersiz |
| 91 | Geçersiz | İsteğe bağlı vertex ID metadata | Geçersiz |
| 210,220,230 | Yalnız (0,0,1) | Yalnız (0,0,1) | Yalnız (0,0,1) |
| 67,410 | Model/paper seçimi | Aynı | Aynı |

LINE WCS, LWPOLYLINE OCS koordinatıdır. Yalnız varsayılan extrusion kabul edilerek OCS dönüşüm belirsizliği kaldırılır. Kod 39, lineweight 370 ve polyline width değerleri **sac t'si değildir**. [Autodesk LINE](https://help.autodesk.com/cloudhelp/2024/ENU/AutoCAD-DXF/files/GUID-FCEF5726-53AE-4C43-B4EA-C84EB8686A66.htm), [Autodesk LWPOLYLINE](https://help.autodesk.com/cloudhelp/2024/ENU/AutoCAD-DXF/files/GUID-748FC305-F3F2-4F74-825A-61F04D757A50.htm)

POLYLINE curve-fit, spline-fit, 3D, mesh ve polyface bayrakları reddedilir. VERTEX geometrik flag'leri sıfır olmalıdır. 2D vertex z=0 olmalı, elevation başlıktan gelmelidir. VERTEX layer'ı yoksa veya "0" ise ana layer devralınır; aksi halde ana layer ile eşit olmalıdır. [Autodesk POLYLINE](https://help.autodesk.com/cloudhelp/2023/ENU/AutoCAD-DXF/files/GUID-ABF6B778-BE20-4B49-9B58-A94E64CEFFF3.htm), [Autodesk VERTEX](https://help.autodesk.com/cloudhelp/2024/ENU/AutoCAD-DXF/files/GUID-0741E831-599E-4CBF-91E1-8ADBCFD6556D.htm)

Kapalı polyline son noktası ilkini tam tekrar ediyorsa bu tek kapanış kaydı çıkarılır, raporlanır ve kapanış kenarı bir kez üretilir. Yaklaşık tekrar tolerans normalizasyonuna tabidir; sıfır kenar doğurursa hata verir. Ardışık iç tekrarlar sessiz silinmez.

### 3.4 Düzlem, birim ve IEEE-754

İzin verilen düzlem WCS z=z0'dır; z0 varsayılan sıfırdır. En uygun düzlem uydurulmaz.

\[
|z-z_0|\le \text{plane_tolerance}.
\]

Varsayılan plane_tolerance=0. Pozitif değer açık kullanıcı iznidir; atılan z sapması raporlanır. Çıkış (x,y)'dir; uzaysal döndürme veya eğik düzlemi düzleştirme yoktur. (0,0,-1) normal de reddedilir.

source_length_unit ve target_length_unit açıkça verilmelidir. unspecified yalnız iki tarafta da unspecified ise kullanılabilir. Desteklenen kesin oranlar: mm=10^-3 m, cm=10^-2 m, in=0.0254 m, ft=0.3048 m.

\[
x'=kx,\qquad y'=ky,\qquad t'=kt.
\]

Kalınlık eşleme değerleri kaynak birimdedir. node_tolerance ve plane_tolerance hedef birimdedir. plane_z kaynak birimindedir ve k ile dönüştürülür.

$INSUNITS varsa açık kaynak birimiyle karşılaştırılır; çelişki reddedilir. Yoksa kullanıcı birimi esastır. Büyüklükten birim tahmin edilmez.

Decimal token ve kesin birim oranı standart kütüphaneyle işlenerek binary64'e tek nihai yuvarlama yapılır. Sonuç taşarsa OverflowError; tam değer sıfır değilken 0'a yuvarlanırsa FloatingPointError verilir. Kaynak decimal'in her basamağının binary64'te temsil edileceği iddia edilmez.

Uzunlukta hypot, yerel referans ve gerektiğinde mantissa/üs ölçekleme kullanılır. Naive dx²+dy² taşması sonlu bir uzunluğu bozamaz. Nihai çekirdek büyüklüğü temsil edilemiyorsa açık istisna korunur.

### 3.5 Kalınlık

ThicknessMap frozen kayıt alanları:

- by_handle: Mapping[str,float], varsayılan boş.
- by_layer: Mapping[str,float], varsayılan boş.
- use_layer_name: bool, varsayılan True.
- default: float | None, varsayılan None.

Öncelik exact handle > exact layer > THICK_<decimal> layer adı > açık default.

THICK_5.0 ve THICK_2e-3 geçerlidir. Sıfır, negatif, NaN, inf veya artçı metin hatadır. Layer karşılaştırması case-sensitive; handle uppercase hexadecimal ile kanoniktir.

Üst öncelikteki bozuk değer için alt kurala düşülmez. Eşleme yoksa GeometryError verilir. Kullanılmayan config handle/layer anahtarları ValueError üretir; layer filtresi dışında kalmış eşleme de kullanılmamıştır.

Polyline kalınlığı bütün kenarlara uygulanır. Width veya vertex bazlı değişken t yoktur. Bölünen segmentler aynı t'yi taşır.

JSON konfigürasyonu bu dört alanlı kapalı sözlüğün karşılığıdır; fiziksel değerler §2.2 F64, handle/layer anahtarları string olur. API'ye verilmeden önce float değerlerine çözülür. Konfigürasyon DXF metninden çalıştırılamaz.

### 3.6 Deterministik kümeleme ve snapping

Mevcut cluster_nodes mutlak mesafe ve geçişli yakınlık bileşenleri kullanır; değiştirilmez. İçe aktarıcı şu daha katı adaptör politikasını uygular:

1. Tüm uçları hedef birime çevir.
2. Koordinatlara göre deterministik sıralama yap; tam eşitlikte signed-zero bitleri sabit tie-break olsun. Handle veya satır sırası temsilci belirlemesin.
3. Sıralı uçlarda cluster_nodes(...,tol=node_tolerance) çağır.
4. Küme çapını denetle. A–B ve B–C yakın, A–C tolerans dışı ise geçişli zinciri belirsiz küme olarak reddet.
5. Temsilci deterministik ilk fiziksel noktadır; ortalama ile yeni nokta icat edilmez.
6. İçe aktarıcı uçları bu temsilciye snap eder. Her özgün uçtan nihai uca hareket toleransı aşamaz.
7. Tüm değişiklikleri raporla; sıfır uzunluk, duplicate veya örtüşme oluşursa reddet.
8. Çekirdeğin kendi cluster_nodes davranışını, JSON okuyucusunu veya eski geometriyi değiştirme.

node_tolerance sonlu, >=0, varsayılan 1e-9 hedef uzunluk birimidir. max(1,scale) gibi birim bağımlı taban eklenmez. Eşdeğer birim dönüşümünde tolerans da k ile ölçeklenmelidir.

Snapping fiziksel geometri değişikliğidir. Değişmezlik normalize geometri için sınanır; ham geometrinin rijitliğiyle koşulsuz bit eşitliği vaat edilmez.

### 3.7 T-birleşimi ve kesişim

junction_policy varsayılan "split_endpoint"; diğer seçenek "reject".

p ucu ve a→b segmenti için analitik izdüşüm:

\[
L=\|b-a\|,\quad \hat t=(b-a)/L,\quad
\ell=(p-a)\cdot\hat t,\quad p_\perp=a+\ell\hat t.
\]

T adayı: 0<ell<L ve uzaklık <=node_tolerance. Uç yakınındaki izdüşüm yeniden bölünmez, uç kümesine tabidir.

- Tam T'de taşıyıcı segment ikiye bölünür.
- Yakın T'de düğüm ve bütün incident uçları izdüşüme taşınır; özgün uca toplam hareket sınırı denetlenir.
- Adaylar aynı geometri anlık görüntüsünden topluca bulunur. Sıraya bağlı ardışık onarım yapılmaz.
- Bir düğüm birden fazla farklı izdüşüm gerektiriyorsa TopologyError; “en yakını seç” yoktur.
- Aynı taşıyıcıdaki bölmeler parametreye göre sıralanır. Eşdeğer bölmeler birleştirilir; sıfır segment yaratılmaz.
- Normalizasyon yeni belirsiz temas yaratırsa sonsuz tekrar yerine hata verilir.
- İç–iç X kesişimleri bu fazda reddedilir; kullanıcı DXF'i ortak uçlarla bölmelidir.
- Pozitif uzunluklu kollinear örtüşme ve ters duplicate daima hatadır.
- reject modunda uç–iç temas da hatadır.

Tam bölmede L=L1+L2; sabit t için alan ve moment integralleri yuvarlanma düzeyinde korunur. Alt segment provenance'ı ebeveyn entity'yi gösterir.

### 3.8 Sınıf seçimi ve hata

Normalize düğümler koordinat sırasına, yönsüz segmentler uç çiftine göre sıralanır. Segment yönü küçük uçtan büyük ucadır. ID'ler n0,n1 ve e0,e1 biçimindedir. Dosya sırası veya kaynak yönü fiziksel ID düzenini değiştirmez.

section_kind:

- open: Section doğrulamasını geçmelidir.
- closed: ClosedSection doğrulamasını geçmelidir.
- mixed: MixedSection olarak kurulur.
- auto: bağlı ağaç → Section; çevrimli ve açık köprüsüz uygun grafik → ClosedSection; çevrim+açık köprü → MixedSection.

Bağlantısız grafik ve boş seçim hatadır. Mevcut topoloji doğrulaması kullanılır; yeni H veya mekanik Tarjan uygulaması yazılmaz.

Sözdizimi/config ValueError; düzlem/yay/uzunluk/t GeometryError; bağlantı/kesişim/örtüşme TopologyError; aritmetik aralık OverflowError/FloatingPointError'dır. Daha özel çekirdek hataları maskelenmez.

### 3.9 API, kaynak sınırı ve provenance

~~~text
loads_dxf(text: str, *, options: DxfImportOptions) -> DxfImportResult
read_dxf(path: str | Path, *, options: DxfImportOptions, encoding: str = "ascii") -> DxfImportResult
~~~

DxfImportOptions tam alanları:

| Alan | Varsayılan/kural |
|---|---|
| source_length_unit | Zorunlu |
| target_length_unit | Zorunlu |
| thickness | Zorunlu ThicknessMap |
| section_kind | auto |
| layers | frozenset[str] veya None; None tüm model-space |
| node_tolerance | 1e-9 |
| plane_z | 0 |
| plane_tolerance | 0 |
| junction_policy | split_endpoint |
| safety_factor | 1e4; çekirdek sözleşmesi |
| max_bytes | 16777216 |
| max_entities | 10000 |
| max_output_segments | 5000 |
| max_pair_checks | 2000000 |

DxfImportResult alanları section ve report'tur. Rapor şunları eksiksiz taşır:

- Sürüm, kaynak/hedef birim, kesin dönüşüm oranı.
- Seçilen/dışlanan entity türleri ve sayıları.
- Entity handle/layer/kenar → çıkış segment provenance'ı.
- İlk/son koordinatlarla snapping ve T-splitting kayıtları.
- En büyük düzlem dışı atım ve düzlem içi hareket.
- Kapanış tekrar kaydı kaldırma sayısı.
- Son düğüm/kenar, bileşen ve çevrim rankı; sınıf.
- Her entity'nin kalınlık kaynağı.

Layer filtresi açık kullanıcı seçimi olup dışlama raporlanır. Paper-space dışlanır ve sayılır; model/paper işaretleri çelişirse hata verilir. Seçim dışı geometry çözümlenmese de dosya yapısı doğrulanır.

ASCII varsayılandır. Legacy layer adları için açık Python text codec seçilebilir; replacement decoding ve encoding tahmini yoktur. $DWGCODEPAGE ile açık codec çelişkisi reddedilir.

O(n²) çift taramalarının maliyeti belgelenir; max_pair_checks aşılmadan önce işlem reddedilir. Kısmi kesit dönülmez. Rapor satır numaralarının sıra bağımlı olması fiziksel determinizm ihlali değildir.

## 4. Teknik Çizim ve Görselleştirme Şartnamesi

### 4.1 Mimari ve API

Matplotlib tembel yüklenir. Yalnız thinwallx veya serialization/dxf import etmek matplotlib.pyplot, bir GUI toolkit veya pencere başlatamaz. Matplotlib kurulu değilken mekanik ve I/O işlevleri kullanılabilir; çizim çağrısı açıklayıcı ImportError üretir.

Bu fazın zorunlu public API'si dosya odaklıdır. Her çağrı kendi Figure'ünü oluşturur, kaydeder ve serbest bırakır; kullanıcıya canlı GUI nesnesi verilmez.

~~~text
plot_geometry(section: Section | ClosedSection | MixedSection, path: str | Path, *, options: GeometryPlotOptions | None = None, units: UnitSystem | None = None, overwrite: bool = False) -> PlotReport
plot_shear_flow(section: Section | ClosedSection | MixedSection, flow: ShearFlowResult | ClosedShearFlowResult | MixedShearFlowResult, path: str | Path, *, options: ShearPlotOptions | None = None, units: UnitSystem | None = None, overwrite: bool = False) -> PlotReport
plot_stresses(section: Section | ClosedSection | MixedSection, result: StressRecoveryResult, path: str | Path, *, options: StressPlotOptions | None = None, units: UnitSystem | None = None, overwrite: bool = False) -> PlotReport
~~~

İmzadaki flow union'ı mevcut modüllerin public sonuç tiplerini ifade eder; yeni mekanik sonuç sınıfı üretilmez. Fonksiyonlar section veya sonuç nesnesini değiştirmez. Aynı kesit üzerinde ardışık çizimler cache dışı fiziksel mutasyon yapamaz.

PlotReport frozen kayıt alanları:

| Alan | İçerik |
|---|---|
| path | Oluşturulmuş mutlak Path |
| format | png veya svg |
| pixel_size | PNG için iki int; SVG için null |
| segment_count | Çizilmiş segment sayısı |
| sample_count | Alan değerlendirme noktalarının toplam sayısı |
| display_origin | İki float; yalnız görüntüleme ötelemesi |
| display_length_exponent | İkili uzunluk ölçek üssü |
| field_scale | Alan yoksa null; varsa ölçek gösterimi |
| skipped_annotations | İstenmemiş/not-applicable açıklamalarının tuple'ı |
| renderer_versions | Python, matplotlib ve numpy sürüm metinleri |

Rapor sayısal mekanik çıktı yerine geçmez. Çizici yeni bir stres maksimumu veya load_factor raporlamaz.

### 4.2 Grafik seçenekleri

Bütün seçenek dataclass'ları frozen ve açık doğrulamalıdır. Bool/int ayrımı, pozitif sonlu boyutlar ve kaynak sınırları denetlenir.

Ortak PlotStyle:

| Alan | Varsayılan |
|---|---|
| figsize | (8.0,6.0) inch |
| dpi | 150, pozitif int |
| font_family | DejaVu Sans |
| font_size | 9.0 point |
| background | white |
| line_width | 1.2 point |
| grid | False |
| max_samples | 200000 |
| svg_hashsalt | thinwallx-v0.8 |
| title | None veya string |

GeometryPlotOptions: style=PlotStyle(), show_thickness=True, show_node_ids=False, show_segment_ids=False, show_centroid=True, show_shear_center=False, show_principal_axes=False, thickness_display_factor=1.0.

ShearPlotOptions: style=PlotStyle(), samples_per_segment=33, arrows_per_segment=3, show_colorbar=True, show_centroid=True, arrow_mode="normalized".

StressPlotOptions: style=PlotStyle(), quantity="sigma_vm", samples_per_segment=65, show_peak=True, abscissa="segment", segment_order=None. quantity seçenekleri sigma_zz, tau_membrane, tau_surface, sigma_vm; abscissa seçenekleri segment ve concatenated.

segment_order verilirse tüm özgün segment indekslerini tam bir kez içeren, değişken uzunluklu integer tuple olmalıdır. Sıralama uzunluğu kesidin segment sayısıyla belirlenir.

Varsayılan stil global rcParams değerlerine bırakılmaz. Tüm etkili renk, font, DPI, line cap/join, arka plan ve SVG ayarları yerel bağlamda belirlenir. Kullanıcının rcParams'ı çağrıdan sonra değişmemiş olmalıdır.

### 4.3 Geometri ve kalınlık bantları

Orta hat:

\[
\mathbf r(\xi)=(1-\xi)\mathbf p_1+\xi\mathbf p_2,\qquad 0\le\xi\le1.
\]

Uçlarda doğrudan özgün uç kullanılır; büyük koordinatta gereksiz yeniden oluşturma yoktur. Yerel koordinatlar üzerinden değerlendirme tercih edilir.

Düz segmentin görsel kalınlık bandı için:

\[
\hat{\mathbf t}=(t_x,t_y),\quad
\hat{\mathbf n}=(-t_y,t_x),\quad
\mathbf r_\pm=\mathbf r\pm \frac{t}{2}\hat{\mathbf n}.
\]

Bandın dört köşesi bu iki uç ofsetinden oluşur. Banda göre alan veya atalet yeniden hesaplanmaz. Birleşim bindirmeleri ve köşe bölgeleri ince cidarlı orta-hat modelinin görsel temsili olup fiziksel katı kesit iddiası taşımaz.

thickness_display_factor pozitif ve sonlu olmalıdır. 1 dışında kullanılırsa grafikte “kalınlık görsel olarak X kat büyütülmüştür” yazısı zorunludur. Gerçek t metadata'da korunur. İnce bandın görünmemesi t'nin artırılması veya mekanik t değişikliği için gerekçe değildir.

- x yatay, y yukarıdır.
- Kesit panelinde equal aspect zorunludur.
- Merkez çizgi koyu gri, bant açık gri ve sabit alpha ile çizilir.
- ID etiketi özgün kullanıcı ID'sini ve gerektiğinde segment indeksini ayırır. Aynı metne sahip ID'ler farklı nesneler sanılmamalıdır.
- Geometrik çizim varsayılan olarak kayma merkezi veya Cw hesaplamaz.
- show_shear_center=True ise mevcut API çağrılır; başarısız hesap sessizce centroid'e çevrilmez.
- İstenmeyen ek bilgiler için tembel hesaplar tetiklenmez.

### 4.4 Centroid, kayma merkezi ve asal eksenler

C siyah dolu daire, S mor çarpı ile işaretlenir; lejant anlamlarını ve birimleri belirtir. Çakışırlarsa iki sembol katmanı kullanılır; ayrı nokta gibi kaydırılmaz.

Mevcut işaret sözleşmesi:

\[
\mathbf I=
\begin{bmatrix}
I_x&-I_{xy}\\
-I_{xy}&I_y
\end{bmatrix},\quad
\theta_p=\frac12\operatorname{atan2}(-2I_{xy},I_x-I_y).
\]

Çizilen birinci doğrultu \((\cos\theta_p,\sin\theta_p)\), ikinci doğrultu \((-\sin\theta_p,\cos\theta_p)\)'dir; çizimde bunların mevcut tensörün I1/I2 özvektörleriyle eşleşmesi test edilir. Aralık veya işaret farklı bir yeni principal-angle formülüyle değiştirilmez.

İzotropik atalette asal doğrultu belirli değildir. Çekirdeğin is_degenerate işareti bu durumu gösterdiğinde keyfi anlamlı bir yön çizilmez; “asal yön izotropi nedeniyle belirsiz” açıklaması konur. Bu işaret otomatik olarak rank-eksik geometri anlamına gelmez.

Asal eksen çizgi uzunluğu görsel bir ölçektir; bir atalet büyüklüğü değildir. C/S görünür sınırın dışında ise panel onları kapsayacak şekilde büyür; gizli kırpma yoktur.

### 4.5 Kayma akışı okları

Akış sonuçları özgün segment sırası, yönü, uzunluğu ve varsa sanal kesim noktalarıyla eşleştirilir. flow.section ile verilen section'ın ham segment koordinatları, t'leri, sırası ve yönü eşleşmelidir; farklı bir kesidin alanı yakıştırılamaz. q(s) mevcut public evaluator üzerinden okunur.

Fiziksel vektör:

\[
\mathbf q(s)=q(s)\hat{\mathbf t}.
\]

Segment ters tanımında scalar q işareti ve koordinatı uygun dönüşür; \(\mathbf q\) aynı fiziksel konumda değişmemelidir. Çizici sadece abs(q) kullanarak yön bilgisini yok edemez.

Ok konumları segmentin içindedir; arrows_per_segment=k için xi_j=(j+1)/(k+1), 0<=j<k integer indeksleridir. Uçlardaki birleşim kalabalığı bu şekilde azaltılır. Sanal kesim noktasında taraf seçimi açık olmalı, fiziksel toplam akış kullanılır; yalnız q_b çizilip toplam q etiketi yazılamaz.

arrow_mode="normalized": ok boyu \(|q|/q_{\max}\) ile sınırlı görsel uzunluk arasında ölçeklenir. q_max=0 ise ok çizilmez ve sıfır alan açıklaması yazılır. Bu gösterim mutlak fiziksel ok uzunluğu değildir; lejant ölçeği belirtir.

Signed q renkleri sıfır merkezli coolwarm paletiyle gösterilir. Aralık [-q_max,+q_max]'tır. Renk çubuğu F/L birimini veya unspecified ifadesini taşır. q_max çizilen alanın analitik uç/ekstremum değerlerini kapsamalıdır; yalnız seyrek renk örneklerine göre kritik değer kırpılamaz.

Sanal kesimli parçalı polinomlar kendi aralıklarında değerlendirilir. Kesimdeki küçük yuvarlama farkı çizicide ortalama ile gizlenmez; fiziksel süreklilik çekirdeğin sonucu olarak doğrulanır.

### 4.6 Gerilme konturu ve 1D dağılım

plot_stresses iki panel üretir: equal-aspect kesit konturu ve 1D segment dağılımı.

Mevcut v0.7 değerlendirmeleri kullanılır:

\[
\sigma_{zz}=a s+b,\quad
\tau_m=c s^2+d s+e,\quad
\tau_{\rm surface}=|\tau_m|+\tau_{\rm sv,surface},
\]
\[
\sigma_{\rm vm}=\sqrt{\sigma_{zz}^2+3\tau_{\rm surface}^2}.
\]

Çizici bu son formülü naive kare toplamıyla yeniden uygulamaz; ölçek korumalı mevcut evaluator'ı kullanır. Katsayı fit'i, yeni gerilme çözümü veya sampled maksimum yoktur.

- sigma_zz ve tau_membrane: sıfır merkezli coolwarm; pozitif/negatif işaret korunur.
- tau_surface ve sigma_vm: sıfırdan başlayan viridis.
- Kontur orta-hat boyunca lif değeridir; gerçek kalınlık boyunca çözümlenmiş 2D gerilme alanı değildir. Bu ifade grafikte veya caption'da bulunur.
- Çizim örnekleri analitik profile'dan alınır. Uçlar, varsa parça sınırları ve kaydedilmiş s_peak mutlaka örnek kümesine dahil edilir.
- Global kritik değer ve konum yalnız StressRecoveryResult'tan alınır.
- Örnek sayısını değiştirmek raporlanan maksimumu değiştiremez.
- Renk normalizasyonu küçük/büyük değerleri güvenli ölçekleyerek yapılır. Bütün değerler sıfırsa sabit sıfır rengi ve “0” renk çubuğu kullanılır; fiziksel 1e-12 tabanı icat edilmez.

### 4.7 Dallanmış grafikte çevre koordinatı

Karma veya dallanmış grafikte bütün segmentleri fiziksel olarak kesintisiz tek s koordinatına sıralamak genel olarak mümkün değildir. Hücre duvarları birden fazla çevrime de ait olabilir.

Varsayılan abscissa="segment": her segment ayrı eğri olarak yerel s∈[0,L_e] üzerinde çizilir; segment ID lejantı bulunur. Aynı x aralığında eğrilerin üst üste gelmesi fiziksel bağlantı iddiası değildir.

abscissa="concatenated": segment_order veya özgün segment sırasıyla yalnız grafiksel koordinat tanımlanır:

\[
s_{\rm plot}=\sum_{j<e}L_j+s_e.
\]

- Eksen etiketi “segment-birleştirilmiş grafik koordinatı; fiziksel sürekli çevre değildir” olmalıdır.
- Her segment sınırında NaN break veya ayrı artist kullanılır; bitiş ile sonraki başlangıç arasında sahte çizgi çizilmez.
- Segment sınırları ve ID'leri işaretlenir.
- Ortak hücre duvarı giriş segmenti olarak bir kez çizilir, her hücre için tekrar eklenmez.
- 1D uzunluk toplamı taşabilecekse görüntüleme ölçeğinde mantissa/üs ile biriktirilir; fiziksel toplam inf olarak sunulmaz.
- Bu fazda Euler yolu bulma, kontur seçme veya fiziksel çevre boyunca yeniden yönlendirme API'si yoktur.

### 4.8 Aşırı ölçeklerde görüntüleme

Çizimde güvenli yerel dönüşüm:

\[
x_{\rm display}=(x-x_{\rm ref})/2^k,\qquad
y_{\rm display}=(y-y_{\rm ref})/2^k.
\]

ref deterministik bir kesit ucudur; k çizilebilir açıklığı normal aralığa getirir. Ara fark taşması da mantissa/üs veya güvenli yerel aritmetikle önlenir. Original değerler ve birimler etiket/raporda korunur. Aynı dönüşüm C, S, bant ve ok konumlarına uygulanır.

- x≈10^12 ötelemede küçük farklar, temsil edilmiş girdide mevcut olduğu ölçüde korunur.
- Girdi yuvarlamasında kaybolmuş koordinat ayrıntısı yeniden üretilemez; çizim bunu “makine hassasiyetinde korunmuş” diye sunamaz.
- Farklı ölçekli segment ekranda piksel altına düşebilir. Bu fiziksel sıfır değildir; raporda küçük ayrıntının çözülemediği belirtilir.
- Alan normalizasyonunda q/scale veya sigma/scale güvenli sırayla hesaplanır; q_max-q_min taşması engellenir.
- Sonlu, geçerli mekanik değerlerin ekran koordinatı için taşması mümkünse yerel ölçek uygulanır; hala gösterilemiyorsa açık ValueError/OverflowError verilir.
- NaN/inf konturu, sessiz kırpma veya DBL_MAX doyurma yoktur.

### 4.9 Headless çalışma ve kaynak ömrü

Matplotlib'in non-interactive Agg çizim yolu ve SVG hardcopy çıktısı kullanılır. GUI backend veya show çağrısı yoktur. Kullanıcının daha önce seçtiği backend zorla değiştirilmez. Önerilen uygulama doğrudan Figure ve FigureCanvasAgg kullanımıdır. [Matplotlib backend belgeleri](https://matplotlib.org/stable/users/explain/figure/backends.html)

Her public çağrı Figure'ünü kendisi oluşturur:

1. Girdileri ve hedef yolu doğrular.
2. Yerel stil bağlamında Figure/canvas oluşturur.
3. Artist'leri ekler, çıktıyı atomik tamamlar.
4. finally içinde sahip olduğu kaynakları serbest bırakır.

Pyplot kullanılmışsa plt.close(fig) zorunludur; close("all") yasaktır. Saf OO Figure manager'a kayıtlı değilse clear ve referansların bırakılması uygulanır; pyplot yalnız kapatma için gereksiz yüklenmez. Başarısız save sırasında da aynı temizlik yapılır. Kullanıcının başka Figure'ü kapanamaz.

### 4.10 PNG/SVG determinizmi ve dosya güvenliği

Aynı girdi, seçenekler ve sabitlenmiş Python/numpy/matplotlib/font/backend ortamında çıktı baytları deterministik olmalıdır.

- Format uzantıdan yalnız .png veya .svg olarak seçilir. Başka uzantı veya format çelişkisi reddedilir.
- PNG boyutu figsize×dpi kuralıyla test edilen tam piksel boyutudur; varsayılan bbox_inches="tight" ile değişken boyut üretilmez.
- SVG viewBox ve boyutlar sabittir.
- Tarih, saat, rastgele UUID, çalışma dizini ve ortamdan gelen kullanıcı adı metadata'ya eklenmez.
- SVG hashsalt sabittir; metadata Date kaldırılır. PNG metadata sabittir.
- Kullanıcı ID/title metni güvenli metin olarak işlenir; TeX çalıştırma kapalıdır, dış font veya URL yüklenmez.
- Aynı süreçte iki yazım ve ayrı temiz süreçlerde yazım hash eşitliğiyle sınanır.
- Farklı işletim sistemi/font/render sürümlerinde bayt eşitliği garanti edilmez. Bu durumda artist koordinatları, etiketler, renk ölçekleri ve sayısal veri eşitliği aranır.
- SVG script ve dış kaynak referansı üretilemez.
- overwrite=False varsayılanı korunur. Çıktı yalnız kullanıcının belirttiği hedefe yazılır.
- Kısmi veya başarısız çıktı geçerli rapor olarak dönülemez.

Metadata ve dosya biçimi davranışı seçilen renderer sürümünde sabitlenip test edilir. [Matplotlib savefig](https://matplotlib.org/stable/api/_as_gen/matplotlib.figure.Figure.savefig.html), [Matplotlib configuration](https://matplotlib.org/stable/users/explain/configuration.html)

## 5. İnvariantlar ve Doğrulama Matrisi

### 5.1 Test prensipleri

Testler mekanik sonucu üretim kodundan kopyalayarak oracle oluşturamaz. Beklenen A, I, J ve gerilme değerleri Fraction veya Decimal ile bağımsız türetilir. Decimal bağlamı en az 100 anlamlı basamaktır; karşılaştırmanın son aşamasında float'a dönülür.

Codec testinde aynı kodlayıcı ve çözücünün aynı hatayı yapmasını önlemek için en az bir elle yazılmış geçerli belge, bir elle hesaplanmış F64 bit deseni ve bozuk belgeler bulunmalıdır. Sadece encode→decode testi yeterli değildir.

Varsayılan np.isclose/pytest.approx mutlak toleransı kullanılamaz. Sıfır olmayan referansta atol=0 ve açık rtol; sıfırda boyutu doğru bağımsız ölçek kullanılır. Hata normları:

\[
E_{\rm rel}=\frac{|f-f_{\rm ref}|}{|f_{\rm ref}|},\quad f_{\rm ref}\ne0,
\]
\[
E_0=\frac{|f|}{S_f},\quad f_{\rm ref}=0,\ S_f>0.
\]

S_f testin analitik büyüklüğünden gelir; max(1,scale) yoktur. Gerçek tüm-sıfır codec testi bit eşitliğiyle yapılır. Aritmetik karşılaştırma sıfıra yakın pozitif sonucu yanlışlıkla kabul edemez.

Varsayılan eşikler:

| Kontrol | Eşik |
|---|---|
| JSON fiziksel skaler aktarımı | Bit eşitliği; 0 hata |
| Yapı, ID türü, sıra, sınıf | Tam eşitlik |
| Normal ölçekli analitik A/I/J/gerilme | rtol<=1e-12, atol=0 |
| Büyük öteleme/dönme, temsili korunmuş girdi | rtol<=1e-10, atol=0 |
| Tam segment bölmenin integral korunumu | Boyutlu referansa göre <=1e-12 |
| Fiziksel sıfır artık | Bağımsız fiziksel ölçeğe göre <=1e-10 |
| PNG/SVG aynı sabit ortam | Tam byte/hash eşitliği |

Daha gevşek eşik gerektiren test otomatik gevşetilemez; gerekçesi ve hata bütçesi incelemeye sunulur.

### 5.2 Bağımsız analitik benchmark'lar

#### 5.2.1 Eş kalınlıklı dikdörtgen orta-hat kutusu

Orta-hat genişliği b, yüksekliği h, kalınlığı t; koordinatlar merkezi orijinde:

\[
A=2t(b+h),\quad C=(0,0),\quad I_{xy}=0,
\]
\[
I_x=t h^2(b/2+h/6),\qquad
I_y=t b^2(h/2+b/6),
\]
\[
A_{\rm cell}=bh,\qquad
J_{\rm BB}=\frac{4b^2h^2}{2(b+h)/t}
=\frac{2tb^2h^2}{b+h}.
\]

Bu ifadeler dört düz segmentin kapalı integrali ve tek hücre Bredt–Batho uyumluluğundan türetilir; üretim sec.Ix veya sec.J oracle içinde çağrılmaz.

N, Mx ve My için:

\[
\sigma(x,y)=N/A+M_y x/I_y-M_x y/I_x.
\]

Tsv için kapalı membran akışı:

\[
q_T=T_{\rm sv}/(2bh),\qquad \tau_T=q_T/t,
\]
\[
\sigma_{\rm vm}=\sqrt{\sigma^2+3\tau_T^2}.
\]

İşaret çevrim yönüne göre eşlenir. Bu benchmark en az bir R12 POLYLINE ve bir AC1015 LINE/LWPOLYLINE dosyasından üretilir. Kare özel durumu J=t b³ olur.

#### 5.2.2 Açık iki bacaklı kesit

(0,0)→(b,0) ve (0,0)→(0,h), aynı t:

\[
A=t(b+h),\quad
x_c=\frac{b^2}{2(b+h)},\quad
y_c=\frac{h^2}{2(b+h)},
\]
\[
I_x=th^3/3-Ay_c^2,\quad
I_y=tb^3/3-Ax_c^2,\quad
I_{xy}=-Ax_cy_c.
\]

Eksenel yükte sigma=N/A sabittir. Genel eğilmede:

\[
D=I_xI_y-I_{xy}^2,
\]
\[
\sigma=N/A+
\frac{M_yI_x+M_xI_{xy}}{D}(x-x_c)
-\frac{M_xI_y+M_yI_{xy}}{D}(y-y_c).
\]

Bütün referanslar Fraction ile hesaplanır. Bu test Ixy işaretini, endpoint yönünü ve centroid çevirisini birlikte denetler. Ortak uçtaki sac bindirmesi analitik modele eklenmez.

#### 5.2.3 T-splitting ve karma topoloji

Bir düz taşıyıcı segmentte tam T noktasıyla bölme, integralin toplamsallığıyla doğrulanır:

\[
\int_0^L f(s)t\,ds
=
\int_0^{L_1}f(s)t\,ds+
\int_{L_1}^{L}f(s)t\,ds.
\]

f=1,x,y,x²,y²,xy için beklenen rasyonel değerler ayrı hesaplanır.

Karma fixture: dikdörtgen kapalı çekirdek ve uçlara bağlı düz açık kollar. Kapalı duvarlar J_open'a dahil edilmez:

\[
J_{\rm total}=J_{\rm BB}+\frac13\sum_{\text{yalnız açık kollar}} L_e t_e^3.
\]

Bu benchmark importer'ın topoloji seçimini ve kalınlık eşlemesini sınar; mevcut karma çözücünün yerine yeni çözüm yazılmaz. Barbell fixture'da iki hücre katkısı ve açık köprü katkısı bağımsız toplanır.

#### 5.2.4 Enine kayma ve süreklilik

§5.2.2 açık L kesidinde iki segment de birleşimden serbest uca yönlendirilsin. Analitik ataletlerden:

\[
\alpha_x=(I_xV_x-I_{xy}V_y)/D,\qquad
\alpha_y=(I_yV_y-I_{xy}V_x)/D.
\]

Yatay bacakta 0<=s<=b için serbest uç tarafının momentleri:

\[
Q_y=t[(b^2-s^2)/2-x_c(b-s)],\qquad
Q_x=-t y_c(b-s).
\]

Düşey bacakta 0<=s<=h:

\[
Q_y=-t x_c(h-s),\qquad
Q_x=t[(h^2-s^2)/2-y_c(h-s)].
\]

Her bacakta q=alpha_x Q_y+alpha_y Q_x. Formül, serbest uç tarafında t integral r_c ds alınmasından gelir; üretim shear-flow fonksiyonu kullanılmaz.

b=h=1, t=1/100, Vx=0, Vy=1 için tam rasyonel örnek:

\[
q_{\rm yatay}(s)=-3/4+3s-9s^2/4,\qquad
q_{\rm düşey}(s)=3/4+3s-15s^2/4.
\]

Her serbest uçta q(1)=0, birleşimde iki dışa yönlü q toplamı sıfırdır. Entegre yatay kuvvet sıfır, düşey kuvvet birdir. Importer canonical yönü farklıysa analitik vektör ilgili yöne dönüştürülür.

Kapalı dikdörtgende alt duvarın soldan sağa yönündeki xi_c konumu sanal kesim seçilsin; diğer üç duvar spanning tree'dedir ve pozitif hücre dolaşımı saat yönünün tersidir. Vy yükünde:

\[
q_0=-\frac{V_y\,t\,b\,h(1-2\xi_c)}{4I_x},
\qquad I_x=t h^2(b/2+h/6).
\]

Bu ifade q0=-integral(q_b/t ds)/integral(1/t ds) koşulundan, dört düz duvarın kesimden başlayan statik momentlerinin kapalı integralleriyle elde edilir. xi_c değerleri 1/6,1/4,9/20,1/2 olarak Fraction kullanılır. Yalnız q0 değil, toplam q'nun kesim iki tarafında sürekliliği ve kuvvet bileşkesi de denetlenir. Canonical segment indeksleri değişebileceğinden kesilen alt duvar koordinatlarıyla seçilir.

Test docstring'i kullanılan yön, kesim ve integral adımlarını gösterir. Testin beklenen değerini üretim kesit ataletinden veya np.linalg.solve'dan üretmek yasaktır.

### 5.3 Öteleme, dönme, sıra ve yön

- Öteleme testi en az 10^12 ve ayrıca ikili üslerle tam temsil edilen örnek içerir.
- Kaynak decimal koordinatlar binary64'e dönüşürken özellik kayboluyorsa test yuvarlanmış gerçek girdiye göre oracle üretir; kaybolmuş geometriyi doğru kabul eden tolerans kullanılmaz.
- Dönmede yük vektörleri de aynı fiziksel dönüşüme tabi tutulur. Sabit sayısal Mx/My ile geometriyi döndürüp aynı gerilmeyi beklemek geçersiz testtir.
- İçe aktarılan kesitlerin farklı entity sırası ve ters yönleri canonical geometride aynı olmalıdır.
- Scalar q işareti değil, aynı noktadaki q t_hat vektörü karşılaştırılır.
- omega* veya shear-center yeniden hesaplanması gereken testlerde mevcut çekirdek API'si kullanılır; adaptör çıktısı düzeltilemez.
- Grafik panelinin baytları öteleme altında aynı olmak zorunda değildir; eksen/origin etiketleri değişebilir. Fiziksel alan ve normalize artist konumları denetlenir.

### 5.4 Test matrisi T01–T60

Her satır en az bir ayrı test veya açık parametrizasyon grubu olarak izlenebilir olmalıdır.

| ID | Test | Zorunlu doğrulama |
|---|---|---|
| T01 | Section JSON round-trip | Koordinat, t, tolerance, sıra, yön, ID türü bit/tam eşit |
| T02 | ClosedSection round-trip | Sınıf, safety_factor, topoloji korunur |
| T03 | MixedSection ve saf-açık MixedSection | Topolojiye bakarak sınıf değiştirilmez |
| T04 | AppliedLoads tüm alanları | İşaret, signed zero, sigma_yield/null korunur |
| T05 | StressRecoveryResult tüm alanları | Katsayılar, profile sırası, peaks ve bileşkeler korunur |
| T06 | F64 uç değerler | ±0, minsub, maxfloat ve yakın komşu float bit eşit |
| T07 | PortableId | int/string/bool/tuple ayrımı; Python key çakışması reddi |
| T08 | Duplicate/unknown/missing JSON alanı | Her nesting düzeyinde açık hata |
| T09 | Bozuk F64 ve NaN/Infinity | Kanonik olmayan metin ve yasak etiket reddi |
| T10 | Geçersiz kesit JSON | Negatif t, ref hatası, sıfır kenar, disconnect reddi |
| T11 | Validate bypass | validate=False kaynak veya wire validate alanı güvenlik açmaz |
| T12 | Birim zarfı | Decode→encode metadata kaybı ve sessiz relabel yok |
| T13 | Result tahrifatı | Yanlış peak/profile key/aralık ve load_factor etiketi reddi |
| T14 | JSON kaynak limitleri | Boyut, derinlik, string ve dizi sınırları |
| T15 | JSON deterministik dosya | Sabit metin, atomik yazma, overwrite koruması |
| T16 | R12 LINE | Bağımsız açık kesit A/C/I/gerilme oracle |
| T17 | R12 POLYLINE | VERTEX/SEQEND, açık ve kapalı zincir |
| T18 | AC1015 LWPOLYLINE | Vertex sayısı, elevation, closed flag |
| T19 | Sürüm ve binary reddi | Eksik/uygunsuz ACADVER, binary imza |
| T20 | Bozuk çift/durum makinesi | Eksik değer/EOF, orphan vertex, tekrarlı entity alanı |
| T21 | Yay ve eğri reddi | ARC/SPLINE ve en küçük sıfır-olmayan bulge dahil |
| T22 | Düzlem/extrusion | Eğik normal, z uyuşmazlığı, 3D flag; açık z0 kabulü |
| T23 | Kalınlık önceliği | Handle > layer > ad > default; bozuk üst değer fallback yapmaz |
| T24 | Width/39/370 semantiği | t sanılmaz; yasak nonzero width/39 reddi |
| T25 | Birim eşdeğerliği | mm↔m, t ve tolerance dönüşümü, INSUNITS çelişkisi |
| T26 | Entity sıra değişimi | Canonical geometri, ID ve fiziksel sonuç aynı |
| T27 | Endpoint/polyline yön tersliği | Fiziksel tensor/gerilme ve q vektörü aynı |
| T28 | Büyük öteleme | Temsili korunmuş geometride 10^12 ve ikili üs fixture |
| T29 | Çok küçük/büyük sayı | Parsing/dönüşümde inf/sessiz sıfır yok |
| T30 | Küme sınırı | tol altı, tam tol, nextafter ile üstü; tol=0 |
| T31 | Geçişli küme zinciri | Çap>tol belirsizliği, sıra bağımsız reddetme |
| T32 | Duplicate/örtüşme | Ters duplicate ve kısmi kollinear örtüşme reddi |
| T33 | Tam T-splitting | A/statik moment/inertia toplamsallığı ve provenance |
| T34 | Yakın/çoklu T | Hareket bütçesi, ambiguity ve reject modu |
| T35 | X kesişimi | Sessiz köprü/bağlılık uydurulmadan hata |
| T36 | Topoloji sınıfı | Open/closed/mixed/auto; açık ekli closed reddi |
| T37 | Barbell/yıldız | İki hücre, köprü, çoklu açık dal ve bağımsız J oracle |
| T38 | Kapalı tekrar ucu | Kapanış kenarı tek; iç tekrar reddi |
| T39 | Katman/model-space | Dışlamalar raporlu; seçili INSERT hata |
| T40 | Encoding/metadata | Legacy layer, bozuk encoding, duplicate handle |
| T41 | Parser kaynak limitleri | O(n²) bütçe aşımı ve kısmi sonuç dönmemesi |
| T42 | Dikdörtgen kutu oracle | A/I/J_BB ve N/M/T altında lif sigma_vm |
| T43 | Enine kayma oracle | Bağımsız integraller, kesim sürekliliği ve bileşke |
| T44 | Matplotlib yokluğu | Çekirdek/JSON/DXF çalışır, plot açıklayıcı hata verir |
| T45 | Headless çalıştırma | Temiz alt süreç, GUI/display olmadan üç plot |
| T46 | Geometri grafiği | Bant köşeleri, equal aspect, ID, C/S doğru |
| T47 | Asal eksen/isotropi | Atalet işareti ve özvektör uyumu; belirsiz eksen etiketi |
| T48 | q ok yönü | Ters segmentte q t_hat aynı; sıfır q ok üretmez |
| T49 | q sanal kesim | Parçalı evaluator, kesim tarafları ve doğru toplam q |
| T50 | Gerilme konturu | Signed/diverging ve vm/sequential renk, gerçek peak işareti |
| T51 | Dallanmış 1D eksen | Segment kopuklukları; ortak duvar tek, sahte süreklilik yok |
| T52 | Çizim aşırı ölçek | Minsub/maxfloat alan, büyük origin, sonlu artist koordinatları |
| T53 | Örnek sayısı bağımsızlığı | 17/65/257 örnekte mekanik maksimum değişmez |
| T54 | PNG/SVG determinizmi | Aynı süreç ve temiz süreçlerde sabit ortam hash eşit |
| T55 | Format doğrulaması | PNG signature/boyut; SVG XML/viewBox/dış kaynak yok |
| T56 | Kaynak ömrü | En az 100 çizim, manager artışı yok; hata yolunda cleanup |
| T57 | rcParams ve nesne mutasyonu | Önce/sonra stil ve kesit/result verisi aynı |
| T58 | Güvensiz title/ID metni | XML escape, TeX/komut/dış dosya yürütme yok |
| T59 | Dosya hata yolları | overwrite=False, izin hatası, yarım çıktı bırakmama |
| T60 | Tüm eski regresyonlar | Tam suite -W error; skip/silme/tolerans gevşetme yok |

### 5.5 Grafik testlerinin niteliği

Sadece “dosya oluştu ve boyutu sıfırdan büyük” testi yeterli değildir.

- PNG, standart kütüphaneyle signature ve IHDR boyutları bakımından denetlenebilir.
- SVG, xml.etree.ElementTree ile parse edilir; boyut, etiket ve tehlikeli dış referans yokluğu denetlenir.
- Artist koordinatları ve renk normalizasyonları küçük iç çizim yardımcılarının testleriyle doğrulanır.
- Aynı ortamda hash testi determinizmi kanıtlar; tek başına fiziksel doğruluğu kanıtlamaz.
- En az açık L, kapalı kutu ve karma barbell çıktıları insan tarafından gözle incelenir; sayısal oracle testleri bu incelemenin yerine bırakılmaz.
- Bellek testi yalnız tek bir RSS değerine bağlanmaz; figure manager sayısı ve tekrarlı çağrıların kalıcı kaynak tutmaması denetlenir.
- Başarısız save ve başarısız mekanik annotation yolu da finally temizliği bakımından test edilir.

## 6. Kabul Kriterleri

### 6.1 Kapsam ve mimari kontrol listesi

- [ ] Yeni işler yalnız v0.8 I/O ve statik çizim kapsamındadır.
- [ ] ACTIVE_PHASE.md kullanıcı talimatı olmadan değiştirilmemiştir.
- [ ] Dondurulmuş mekanik kaynaklar, formüller, toleranslar ve test beklentileri değişmemiştir.
- [ ] Sıfır CAD/simgesel/FEM bağımlılığı korunmuştur.
- [ ] Matplotlib yalnız isteğe bağlı çizim yolunda yüklenmektedir.
- [ ] Python >=3.10 ve açık tip ipuçları sağlanmıştır.
- [ ] Yeni nontrivial formül/algoritmalar kaynak veya türetim notu taşımaktadır.
- [ ] İçe aktarma ve çizim fiziksel sonuçları sessiz düzeltmemektedir.

### 6.2 JSON kontrol listesi

- [ ] Beş hedef nesne tipi için bütün tanımlanmış kamusal veriler round-trip olmaktadır.
- [ ] Fiziksel float'lar normal/alt-normal/signed-zero dahil bit eşittir.
- [ ] JSON Schema ve Python doğrulayıcı aynı kapalı alan sözleşmesini uygulamaktadır.
- [ ] Unknown/duplicate/missing alan, belirsiz ID ve tanımsız ref reddedilmektedir.
- [ ] Birim metadata'sı kaybolmamaktadır.
- [ ] validate bypass yolu bulunmamaktadır.
- [ ] NaN/inf sızıntısı yoktur; yalnız semantik yük faktörü etiketi vardır.
- [ ] Sonuç tutarlılığı ile yeniden mekanik doğrulama arasındaki sınır belgelenmiştir.
- [ ] Dosya yazma atomik ve varsayılan olarak overwrite korumalıdır.
- [ ] Boyut/derinlik limitleri güvenilmeyen girdiye uygulanmaktadır.

### 6.3 DXF kontrol listesi

- [ ] Desteklenen sürüm/entity matrisi eksiksiz test edilmiştir.
- [ ] LINE WCS ve polyline OCS/elevation ayrımı doğrudur.
- [ ] Kod 39, width ve lineweight sac kalınlığına dönüştürülmemektedir.
- [ ] Katman/handle kalınlık önceliği ve kaynak birimi açıktır.
- [ ] Düzlem dışı, eğri, mesh ve blok durumları açıkça reddedilmektedir.
- [ ] Cluster representative ve T-splitting sonuçları sıra bağımsızdır.
- [ ] Geçişli tolerans zinciri, duplicate, örtüşme ve X kesişimi test edilmiştir.
- [ ] Bağlantısızlık ve sınıf/topoloji uyumsuzluğu sessiz kabul edilmemektedir.
- [ ] Tüm geometri değişiklikleri provenance raporunda görünmektedir.
- [ ] Büyük/küçük ölçeklerde aritmetik taşma ve sıfırlanma açık hatadır.
- [ ] Rasyonel/100 basamak oracle'lar üretim özelliklerinden bağımsızdır.

### 6.4 Grafik kontrol listesi

- [ ] Üç public plot API PNG ve SVG üretmektedir.
- [ ] GUI, show, event-loop ve dış TeX/font çalışma zamanı yoktur.
- [ ] Equal aspect, kalınlık bantları, C/S ve principal-axis işaretleri doğrudur.
- [ ] q okları fiziksel q t_hat yönünü göstermektedir.
- [ ] Gerilme konturu analitik evaluator kullanmaktadır.
- [ ] Kritik değer sampled grid'den değil v0.7 sonucundan alınmaktadır.
- [ ] Dallanmış kesitlerin 1D ekseni sahte fiziksel çevre sunmamaktadır.
- [ ] Görsel ölçekler ve kalınlık büyütmesi açıkça etiketlenmiştir.
- [ ] Sıfır ve aşırı ölçekli alanlar inf/nan artist üretmemektedir.
- [ ] Sabit ortamda PNG/SVG bayt determinizmi doğrulanmıştır.
- [ ] Başarı ve hata yollarında Figure kaynakları temizlenmektedir.
- [ ] Kullanıcı rcParams'ı ve diğer figürler değişmemektedir.

### 6.5 Zorunlu test çalıştırma ve dondurma kararı

Uygulama tamamlandığında önce hedef testler:

~~~console
python -m pytest -W error tests/test_serialization.py tests/test_dxf_import.py tests/test_dxf_benchmarks.py tests/test_plotting.py
~~~

Sonra filtrelenmemiş tam paket:

~~~console
python -m pytest -W error
~~~

Bu komutlar şartnamenin gelecekteki uygulama kabul prosedürüdür; bu belgenin yazılması sırasında çalıştırıldıkları veya test dosyalarının mevcut olduğu iddia edilmez.

Dondurma raporu şunları içerir:

1. Değişen dosyalar ve frozen kaynakların korunma kontrolü.
2. Gerçek test sayıları, komutlar, sıfır uyarı sonucu ve başarısız/atlanan test gerekçeleri.
3. Bağımsız benchmark hata değerleri; kullanılan atol/rtol.
4. JSON bit eşitliği, DXF normalize geometri ve grafik determinizmi sonuçları.
5. Üç örnek PNG/SVG'nin gözle doğrulaması.
6. Kaynak limitleri, desteklenmeyen DXF profilleri ve çapraz renderer determinizmi sınırı.
7. Açık P0/P1/P2 bulgular ve kalan işler.

Herhangi bir sessiz geometri kaybı, yanlış t eşlemesi, mekanik işaret hatası, round-trip veri kaybı, inf/nan sızıntısı, desteklenmeyen geometrinin kabulü veya dondurulmuş çekirdek regresyonu dondurmayı engeller.

Bütün zorunlu kriterler doğrulanmadan “v0.8 dondurulmaya hazırdır” ifadesi kullanılamaz. Dondurma sonrasında STOP: sonraki faza, yeni mekanik modele veya ek CAD desteğine kendiliğinden geçilmez.
