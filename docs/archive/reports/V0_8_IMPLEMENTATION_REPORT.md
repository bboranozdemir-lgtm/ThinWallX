# ThinWallX v0.8 — Uygulama ve Doğrulama Raporu

## Sonuç

JSON, dar profilli ASCII DXF içe aktarımı ve headless teknik çizim modülleri eklendi. Son filtrelenmemiş test koşusu:

~~~console
python -m pytest -W error
520 passed in 20.09s
~~~

Yeni faz testleri 131 adettir; eski 389 test korunmuştur. Son hedef koşusunda 131 test geçti. Ardından yapılan son doğrulama/stil değişiklikleri de yukarıdaki 520 testlik tam koşuda kapsandı. Başarısız veya atlanan test yoktur.

v0.8, kullanıcının açık ve tam onayıyla resmen FROZEN olarak kaydedildi. Paket sürümü pyproject.toml ve thinwallx.__version__ içinde 0.8.0 olarak eşleştirildi. ACTIVE_PHASE.md ve hesaplama algoritmaları değiştirilmedi. Bu onay bağımsız adversarial denetimin yerine geçmez ve matematiksel hatasızlık iddiası taşımaz. Harici paket yayını yapılmadı; çalışma alanında Git deposu olmadığından Git tag oluşturulmadı.

## Değişen/eklenen dosyalar

- src/thinwallx/serialization.py: katı JSON, F64, PortableId, zarf ve atomik dosya yazımı.
- src/thinwallx/dxf.py: ASCII durum makinesi, sürüm/entity doğrulaması, kalınlık eşleme, snapping/T-splitting ve provenance.
- src/thinwallx/plotting.py: geometri, kayma akışı, gerilme konturu ve ayrık 1D grafik.
- schemas/thinwallx-0.8.schema.json: yapısal şema.
- pyproject.toml: yalnız isteğe bağlı plots/matplotlib bağımlılığı eklendi.
- tests/test_serialization.py: 44 test.
- tests/test_dxf_import.py: 53 test.
- tests/test_dxf_benchmarks.py: 9 test.
- tests/test_plotting.py: 25 test.
- tests/fixtures/dxf: R12 açık L, R12 kapalı dikdörtgen, AC1015 barbell.
- examples/v08_io_and_plots.py: tekrarlanabilir örnek üretici.
- examples/v08_output: 9 JSON ve 18 PNG/SVG çıktı.
- docs/V0_8_USAGE.md: kullanım ve sınırlar.
- docs/V0_8_IMPLEMENTATION_REPORT.md: bu rapor.

Uygulama tesliminde mevcut 44 Python kaynak/test dosyasının SHA-256 özetleri karşılaştırıldı: değişiklik yoktu. Sonraki resmî dondurma işleminde bunlar içinde yalnız src/thinwallx/__init__.py sürüm metadatası 0.8.0 olarak değiştirildi; mekanik algoritmalar ve testler korunmuştur. Sürüm güncellemesinden sonraki tam koşu da yukarıdaki 520 başarılı testle tamamlandı. ACTIVE_PHASE.md SHA-256 değeri:

~~~text
5B6C6210F9D527A373CCC2D791F458A791C0816FE14155B54878005F40A928F6
~~~

## Matematiksel ve sayısal yöntemler

- F64 serileştirmesi float.hex/fromhex ile bit düzeyindedir; signed zero ve minimum subnormal korunur.
- DXF decimal ve birim oranları Fraction üzerinden tek binary64 dönüşümüne gider.
- İzdüşüm, yön ve tolerans kareleri exact Fraction aritmetiğiyle sınanır; geometri onarımı sonuç raporunda görünür.
- Mevcut cluster_nodes ve mevcut topology/section kurucuları yeniden kullanılır; çekirdek mekanik formülü eklenmez.
- Çizimde yerel origin ve ikili ölçek kullanılır.
- q renk aralığı, düz/parçalı quadratic alanın uç ve türev kökü adaylarını kapsar.
- Gerilme maksimumu v0.7 sonucundan alınır; örnekleme yalnız çizim içindir.
- Agg/OO Figure yolu, yerel rc bağlamı, atomik kaydetme ve finally temizliği kullanılır.

## Benchmark ve güvence sonuçları

- Açık L için bağımsız Fraction A/Ixy/eğilme gerilmesi referansları geçti.
- Dikdörtgen için A, Ix, Iy, J_BB ve N/M/T altında sigma_vm referansları rtol=1e-12, atol=0 ile geçti.
- Açık L enine q polinomları ve barbell bağımsız J toplamı geçti.
- Kapalı kutuda 1/6, 1/4, 9/20, 1/2 sanal kesim oracle'ları geçti. Fiziksel akış sürekliliği 1e-10 fiziksel ölçekle sınandı.
- Öteleme, ters yön, sıralama, T-bölme, geçişli tolerans zinciri, örtüşme ve X kesişimi testleri geçti.
- JSON bit testleri +0, -0, 5e-324, maxfloat ve komşu float değerlerini kapsıyor.
- Çizim alan testleri 5e-324/maxfloat; geometri testleri ikili -200/+200 ölçekleri kapsıyor.
- Aynı süreçte ve temiz alt süreçlerde PNG/SVG byte determinizmi geçti.
- 100 tekrarlı çizimde figure-manager artışı görülmedi; save hatasında geçici dosya temizliği geçti.
- Matplotlib import edilmeden çekirdek/JSON/DXF çalışması ve matplotlib yokken açıklayıcı hata doğrulandı.
- Açık L, kutu ve barbell geometri/gerilme çıktıları görsel olarak incelendi. Görülen başlık ve renk çubuğu kırpılması düzeltildi; son kutu kayma çıktısı yeniden incelendi.
- Örnek üretici ayrıca Python uyarıları hata kabul edilerek çalıştırıldı.

## Sınırlar ve kalan yayın işleri

- Test ortamı Windows, Python 3.14.6, pytest 8.4.2, matplotlib 3.11.1'dir. Python 3.10 dahil bütün desteklenen sürüm kombinasyonlarında CI matrisi bu oturumda çalıştırılmadı.
- Genel CAD desteği yoktur: yaylar, eğik OCS, INSERT/blok açılımı, 3D ve variable-width geometri reddedilir.
- ASCII varsayılanıyla legacy codepage çelişkisi açık hatadır; kullanıcı codec seçmelidir.
- Snapping geometri değişikliğidir; ham geometriyle koşulsuz rijitlik eşitliği vaat edilmez.
- Girdi binary64'e dönüşürken kaybolan koordinat ayrıntısı geri üretilemez.
- Renderer ortamı değişirse byte determinizmi garanti edilmez.
- Tek başına stress-result JSON dosyası geometri/yük bağıntısının mekanik sertifikası değildir.
- Kaynaklar içinde çözümsüz TODO/NotImplemented placeholder bırakılmadı.
- Kullanıcı resmî dondurmayı ve 0.8.0 sürümünü onayladı. Bağımsız denetim ve harici yayın bu işlemin kapsamı dışındadır. v0.9 için geçiş izni verilmiştir; mevcut yol haritasında v0.9 kapsamı tanımlı olmadığından uygulamasına başlanmadı.
