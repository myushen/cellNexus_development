#!/usr/bin/env bash
set -euo pipefail
set +x

FILE="/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/hca2025_pseudobulk_se.h5ad"
DEPOSITION="22700163"
ZENODO_ENDPOINT="https://zenodo.org"
MAX_ATTEMPTS=5

: "${ZENODO_TOKEN:?Please set ZENODO_TOKEN}"

FILENAME="$(basename "$FILE")"

BUCKET="$(
  curl --silent --show-error --fail \
    --header "Authorization: Bearer ${ZENODO_TOKEN}" \
    "${ZENODO_ENDPOINT}/api/deposit/depositions/${DEPOSITION}" |
  jq --raw-output '.links.bucket'
)"

for attempt in $(seq 1 "$MAX_ATTEMPTS"); do
  RESPONSE="$(mktemp)"

  echo "Attempt ${attempt}/${MAX_ATTEMPTS}: ${FILENAME}"

  if curl \
    --http1.1 \
    --progress-bar \
    --fail-with-body \
    --connect-timeout 60 \
    --keepalive-time 30 \
    --header "Authorization: Bearer ${ZENODO_TOKEN}" \
    --output "$RESPONSE" \
    --upload-file "$FILE" \
    "${BUCKET}/${FILENAME}"
  then
    echo
    echo "Upload succeeded."
    jq . "$RESPONSE"
    rm -f "$RESPONSE"
    exit 0
  fi

  status=$?
  echo
  echo "Attempt failed, curl status: ${status}"
  cat "$RESPONSE" || true
  rm -f "$RESPONSE"

  # 逐渐延长等待，最高10分钟
  delay=$((attempt * 60))
  (( delay > 600 )) && delay=600

  echo "Retrying from the beginning in ${delay} seconds..."
  sleep "$delay"
done

echo "Upload failed after ${MAX_ATTEMPTS} attempts." >&2
exit 1
