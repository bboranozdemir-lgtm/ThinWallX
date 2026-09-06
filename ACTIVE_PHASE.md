# Sectalix — Active Scope v1.0.1

STATUS: ACTIVE — Sectalix v1.0.1 rename and release preparation only

## Yetkilendirilmiş v1.0.1 istisnası (öncelikli)

Kullanıcı Sectalix → Sectalix yeniden adlandırmasını açıkça onaylamıştır.
İzin yalnızca src/thinwallx → src/sectalix taşıması, import/path/marka/CLI adı,
1.0.1 sürümü, metadata, LICENSE, README, docs, tests, schemas, MANIFEST,
GitHub Actions/OIDC hazırlığı ve bunların test/build/install doğrulamasını kapsar.
Matematiksel/nümerik davranış ve mevcut v0.8 dosya biçimi uyumluluğu korunacaktır.
Commit, push, remote rename/settings, eski tag/release değiştirme, yeni tag ve
gerçek PyPI yayını yasaktır. Önceki push yetkisi bu çalışma için geçersizdir.
Aşağıdaki v1.0 korumaları yalnızca bu dar istisnanın dışında aynen geçerlidir;
tarihsel v1.0.0 etiketi ve algoritmalar dondurulmuş olarak korunur.

Normatif şartname: [docs/archive/specifications/ACTIVE_PHASE_v1.0.md](docs/archive/specifications/ACTIVE_PHASE_v1.0.md). Uygulamadan önce bu belgenin tamamı okunmalıdır.

- Kapsam: Sectalix 1.0.1 yeniden adlandırması, Markdown hesap raporu, teknik dokümantasyon ve paketleme.
- v1.0 kabul kriterleri 559 test ile sağlanmıştır. src/ ve tests/ bütünüyle dondurulmuştur.
- Yeni mekanik formül, GUI veya kapsam dışı bağımlılık eklenemez.
- Önceki aktif v0.7 şartnamesi [docs/archive/specifications/ACTIVE_PHASE_v0.7.md](docs/archive/specifications/ACTIVE_PHASE_v0.7.md) dosyasında korunmuştur.
- Kullanıcı v1.0 dondurma, yerel Git mühürleme ve Yol A dağıtım/CI hazırlığını açıkça onaylamıştır.
- Yol A yalnızca release belgeleri, Git ve CI/CD altyapısını kapsar; mekanik/test değişiklikleri yasaktır.
- Doğrulama: python -m pytest -W error; beklenen 559 passed, sıfır hata ve uyarı.
- Dondurma harici yayın anlamına gelmez. GitHub/PyPI gönderimleri kullanıcıya bırakılır; yeni faz başlatılmaz.

## Yetkilendirilmiş bakım istisnası

Bu bölüm v1.0.0 bakım tarihçesidir; yukarıdaki dar v1.0.1 yeniden adlandırma ve
release-hazırlığı yetkisi önceliklidir ve onun sınırlarını genişletmez.

Kullanıcı mevcut denetim bulgularını gidermek üzere testlerde dar kapsamlı bakım
değişikliklerini ve GitHub main dalına gönderimi onaylamıştır. Bu istisna yalnızca
ölçek-duyarlı tork toleransı, paketleme uyarı filtresi ve ilgili doğrulama/belge/CI
bakımını kapsar. src/ değiştirilemez; yeni mekanik özellik eklenemez. v1.0.0 etiketi
taşınmaz, PyPI yayını yapılmaz ve yeni faz başlatılmaz.
