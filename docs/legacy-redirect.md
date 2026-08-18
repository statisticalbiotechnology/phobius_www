# Redirecting phobius.sbc.su.se to the new service

Notes for whoever administers the old host and the DNS record. Replace
`PHOBIUS-APP` below with the actual Serve subdomain before sending.

## Current state

Verified 2026-08-18:

| | |
|---|---|
| `phobius.sbc.su.se` resolves to | `130.237.65.184` = `es2.scilifelab.se` |
| Server | Apache, serving a static "Phobius - Maintenance" page |
| `/cgi-bin/predict.pl` | returns 404 — the legacy CGI is already gone |
| New service | `https://PHOBIUS-APP.serve.scilifelab.se` |

Two things worth noting before choosing an option:

* The hostname is in the `su.se` zone, so **Stockholm University controls the DNS
  record**, but it already points at **SciLifeLab** infrastructure. A redirect on
  the existing host therefore needs no involvement from SU at all.
* Nothing is currently served from the old address except a maintenance page, so
  either change can be made without an outage.

## Option A — keep the old address working (recommended)

Phobius has been cited with this URL since 2004 and it appears in published
methods sections and in scripts. Keeping it alive preserves all of that.

1. Email <serve@scilifelab.se>, ask them to add `phobius.sbc.su.se` as a custom
   domain for the app, and request the DNS settings. Their documentation says:
   *"you will need to purchase a domain name yourself and set up DNS settings
   that you will get from us."* Here the domain already exists, so only the DNS
   change is needed.
2. Send the values they return to SU IT to apply to the `su.se` zone. This is
   normally a `CNAME`, which means the current `A` record is replaced.
3. Confirm SciLifeLab Serve issues a TLS certificate for `phobius.sbc.su.se`;
   without one the old address will fail with a certificate warning.

Do these in order — the DNS change should not be made before Serve is ready to
answer for the name, or the address will break in the interval.

## Option B — permanent redirect from the existing host

Simpler, and needs nothing from SU. On `es2.scilifelab.se`, in the Apache
configuration for this virtual host:

```apache
<VirtualHost *:443>
    ServerName phobius.sbc.su.se

    # 308, not 301: a permanent redirect preserves the request method and body,
    # so scripts still POSTing to /cgi-bin/predict.pl keep working against the
    # new service, which serves that path for backwards compatibility. A 301 or
    # 302 turns those POSTs into GETs and the submitted sequence is lost.
    RedirectPermanent / https://PHOBIUS-APP.serve.scilifelab.se/
    Redirect 308 / https://PHOBIUS-APP.serve.scilifelab.se/
</VirtualHost>
```

Use **one** of those two lines: `Redirect 308` if old POST clients should keep
working, `RedirectPermanent` (301) if only browser traffic matters. The path is
preserved in both cases, so `/cgi-bin/predict.pl` and `/instructions.html` land
on their equivalents.

The trade-off against Option A is that the address in the browser changes to
`serve.scilifelab.se`, and published references to `phobius.sbc.su.se` then
depend on this redirect staying in place indefinitely.

## Checking afterwards

```bash
curl -sSI https://phobius.sbc.su.se/ | head -3
curl -sS -X POST https://phobius.sbc.su.se/cgi-bin/predict.pl \
     -F 'protseq=>test
MYGKIIFVLLLSAIVSISASSTTGVAMHTSTSSSVTKSYISSQTNDTHKRDTYAATPRAHEVSEISVRTVYPPEEETGERVQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKKSPSDVKPLPSPDTDVPLSSVEIENPETSDQ' \
     -F 'format=short' -L | grep -A2 'SEQENCE ID'
```

Under Option A the first command returns `200`; under Option B it returns `308`
(or `301`) with a `Location:` header. The second should print a prediction line
either way — that is the check that old scripts still work.
