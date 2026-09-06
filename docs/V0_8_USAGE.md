# ThinWallX v0.8 — kullanım

Bu sürüm adaptör modülleri ekler; v0.1–v0.7 mekanik kodu değişmez. API'leri doğrudan ilgili modülden import edin. v0.8 kullanıcı onayıyla resmen dondurulmuştur; paket sürümü 0.8.0'dır.

## JSON

~~~python
from thinwallx import Node, Segment, Section, AppliedLoads
from thinwallx.serialization import (
    UnitSystem, to_dict, from_dict, to_json, from_json, write_json, read_json,
)

section = Section([
    Segment(Node(0.0, 0.0), Node(2.0, 0.0), 0.01, id="flange"),
    Segment(Node(0.0, 0.0), Node(0.0, 1.0), 0.01, id="web"),
])
units = UnitSystem(length="mm", force="N")
wire = to_json(section, units=units)
document = from_json(wire)
restored_section = document.value
assert document.units == units
assert restored_section.segments == section.segments
~~~

from_dict/from_json/read_json bir DecodedDocument döndürür: value ve units. to_dict/to_json/write_json hem ham nesne hem bu zarfı kabul eder. Birim etiketi sayı dönüştürmez; mevcut zarfı farklı birimle yeniden etiketleme reddedilir.

Fiziksel değerler kanonik float.hex string kaydıdır; decimal JSON number kabul edilmez. Signed zero ve alt-normal sayılar bit düzeyinde korunur. JSON şeması schemas/thinwallx-0.8.schema.json dosyasındadır; çapraz referans/topoloji/peak doğrulaması Python okuyucusunda uygulanır.

Desteklenen ID'ler None, str, int, bool, sonlu float ve bu tiplerden tuple'dır. Özel Python nesneleri taşınamaz. Result sözlüğü ID'leri string'e çevrilmez.

AppliedLoads ve StressRecoveryResult aynı codec ile saklanabilir. Sonuç dosyası geometri ve yükleri içermez; bunlar ayrı belgeler olarak saklanmalıdır. Peak katsayı tutarlılığı denetlenir; bağımsız arşiv mekanik sertifika veya dijital imza değildir.

write_json varsayılan olarak mevcut dosyayı ezmez. overwrite=True yalnız açıkça istendiğinde kullanılır. Okuyucular byte/derinlik/kayıt limitleri uygular.

## ASCII DXF

~~~python
from thinwallx.dxf import DxfImportOptions, ThicknessMap, read_dxf

options = DxfImportOptions(
    source_length_unit="mm",
    target_length_unit="mm",
    thickness=ThicknessMap(by_layer={"WEB": 2.0}, default=1.0),
    node_tolerance=1e-9,
    section_kind="auto",
)
imported = read_dxf("section.dxf", options=options)
section = imported.section
report = imported.report
~~~

THICK_5.0 gibi layer adları desteklenir. Öncelik handle > layer > THICK adı > açık default'tur. Eşleme kalınlıkları kaynak biriminde, node_tolerance hedef birimindedir. Kullanılmayan eşlemeler yazım hatasını gizlememek için reddedilir.

AC1009: LINE ve düz 2D POLYLINE. AC1015: ayrıca LWPOLYLINE. Bulge, nonzero width/extrusion thickness, eğik extrusion normal, 3D polyline, INSERT ve yaylar reddedilir. DXF 39/370 kodları sac t'si değildir.

Varsayılan normalize etme endpoint snapping ve açık T-birleşimi bölmesidir. Değişimler report.changes/provenance içindedir. Geçişli yakınlık kümesinin çapı toleransı aşıyorsa hata verilir. X kesişimleri ve örtüşmeler onarılmaz. Sadece exact ham geometri isteyen kullanıcı node_tolerance=0 ve junction_policy="reject" seçebilir.

varsayılan okuma ascii'dir; legacy layer adları için encoding="cp1254" gibi açık codec verilebilir. Codepage çelişkisi, eksik EOF/sürüm ve bozuk grup çiftleri reddedilir. Bu bir genel CAD okuyucusu değildir.

Bağlı açık ağaç Section, saf hücreli grafik ClosedSection, hücre+köprü MixedSection olur. section_kind ile açıkça sınıf istenebilir; uyumsuz geometri hata verir. Grafik işlem bütçesi max_pair_checks ile sınırlıdır; varsayılan 2 milyon muhafazakâr işlem tahminidir.

## Teknik çizimler

Grafik bağımlılığı isteğe bağlıdır:

~~~console
python -m pip install -e ".[plots]"
~~~

~~~python
from thinwallx import AppliedLoads
from thinwallx.plotting import (
    plot_geometry, plot_shear_flow, plot_stresses,
    GeometryPlotOptions, StressPlotOptions,
)

plot_geometry(section, "geometry.svg",
              options=GeometryPlotOptions(show_shear_center=True, show_principal_axes=True))
flow = section.calculate_shear_flow(vx=0.0, vy=1.0)
plot_shear_flow(section, flow, "shear.png")
stress = section.calculate_stresses(AppliedLoads(N=1.0, Mx=0.1))
plot_stresses(section, stress, "stress.svg",
              options=StressPlotOptions(quantity="sigma_vm", abscissa="segment"))
~~~

Çağrılar kendi figürlerini oluşturur ve kapatır; pencere açmaz. Canlı Figure yerine PlotReport döner. Varsayılan overwrite=False'dur. PNG/SVG desteklenir.

Geometri equal-aspect'tir; fiziksel farklar yerel origin ve 2'nin kuvvetiyle normalize edilir. Renk çubuğundaki ölçek görseldir, fiziksel değer değiştirilmez. Tam field_scale raporda float.hex biçimindedir. Kalınlık bantları katı kesit değil, orta-hat modelinin görsel yardımcısıdır.

Kayma okları q*tangent yönündedir. Stress maksimumu v0.7 sonucundan alınır, grafik örnekleriyle aranmaz. abscissa="segment" yerel s eğrileri; "concatenated" ayrık segmentleri grafiksel olarak yan yana koyar ve fiziksel çevre olmadığını etiketler.

Bayt determinizmi aynı Python/numpy/matplotlib/font ortamı için sınanır; farklı render sürümlerinde aynı bayt garantisi yoktur.

## Çalışan örnekler ve testler

examples/v08_io_and_plots.py, tests/fixtures/dxf altındaki açık L, kapalı kutu ve barbell dosyalarını okuyarak 9 JSON ve 18 PNG/SVG üretir. examples/v08_output altında üretilmiş örnekler bulunur. Gerekirse src dizinini PYTHONPATH'e ekleyin veya editable kurulum kullanın.

~~~console
python examples/v08_io_and_plots.py --output examples/v08_output
python -m pytest -W error tests/test_serialization.py tests/test_dxf_import.py tests/test_dxf_benchmarks.py tests/test_plotting.py
python -m pytest -W error
~~~

Hiçbir test sayısı tek başına bütün matematiksel hataların yokluğunun kanıtı değildir. Kullanıcının resmî dondurma onayı bağımsız adversarial denetim yapılmış olduğu iddiasını taşımaz.
