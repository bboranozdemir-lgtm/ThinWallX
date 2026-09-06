# ThinWallX — Active Phase v1.0

STATUS: FROZEN — v1.0 Production-Quality Technical Release Complete

Normatif şartname: [ACTIVE_PHASE_v1.0.md](ACTIVE_PHASE_v1.0.md). Uygulamadan önce bu belgenin tamamı okunmalıdır.

- Kapsam: argparse CLI, Markdown hesap raporu, teknik dokümantasyon ve 1.0.0 paketleme.
- v1.0 kabul kriterleri 559 test ile sağlanmıştır. src/ ve tests/ bütünüyle dondurulmuştur.
- Yeni mekanik formül, GUI veya kapsam dışı bağımlılık eklenemez.
- Önceki aktif v0.7 şartnamesi [ACTIVE_PHASE_v0.7.md](ACTIVE_PHASE_v0.7.md) dosyasında korunmuştur.
- Kullanıcı v1.0 dondurma, yerel Git mühürleme ve Yol A dağıtım/CI hazırlığını açıkça onaylamıştır.
- Yol A yalnızca release belgeleri, Git ve CI/CD altyapısını kapsar; mekanik/test değişiklikleri yasaktır.
- Doğrulama: python -m pytest -W error; beklenen 559 passed, sıfır hata ve uyarı.
- Dondurma harici yayın anlamına gelmez. GitHub/PyPI gönderimleri kullanıcıya bırakılır; yeni faz başlatılmaz.
