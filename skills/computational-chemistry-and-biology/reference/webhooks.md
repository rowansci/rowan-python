# Webhooks

Get notified when a workflow completes instead of polling. Pass `webhook_url=` to any `submit_*_workflow`. When the run finishes, Rowan POSTs the result to that URL with an `X-Rowan-Signature` header (HMAC-SHA256) so your server can verify the request genuinely came from Rowan.

## Manage the signing secret

```python
secret = rowan.create_webhook_secret()   # idempotent, returns the existing secret if already set
secret = rowan.get_webhook_secret()      # current secret, or None if unset
secret = rowan.rotate_webhook_secret()   # new secret, invalidates the old one
```

## Submit with a callback

```python
rowan.create_webhook_secret()            # once, ahead of time
wf = rowan.submit_admet_workflow(
    initial_smiles="...",
    folder=folder,
    webhook_url="https://your-server.com/webhook",
)
```

## Verify an incoming POST (on your server)

When the workflow completes, Rowan POSTs to `webhook_url` with the result body and an `X-Rowan-Signature` header. Verify it before trusting the payload:

```python
ok = rowan.verify_webhook_secret(
    raw_body,             # raw (unparsed) request body bytes
    signature_header,     # value of the X-Rowan-Signature header
    secret,               # your webhook secret
    max_age_seconds=300,  # reject requests older than this (default 5 min)
)
```

Use the **raw** request-body bytes. Re-serializing parsed JSON will change the bytes and break signature verification.
