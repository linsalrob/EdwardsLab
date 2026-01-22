
"""
Search Bluesky for posts containing the hashtag #ABACBS2025
and save the results to a TSV file.

Requires:
    pip install atproto

Set:
export BSKY_HANDLE="your_handle.bsky.social"
export BSKY_APP_PASSWORD="your-app-password-here"


Search Bluesky for posts containing a given hashtag
and save the results to a TSV file.

Works with the current atproto package (2024–2025 API).
"""

import os
import sys
import time
import csv
from atproto import Client

HASHTAG = "ABACBS2025"
OUTPUT_TSV = "bluesky_ABACBS2025.tsv"
MAX_POSTS = 2000
PAGE_LIMIT = 100
SLEEP = 0.3


def get_client():
    handle = os.getenv("BSKY_HANDLE")
    passwd = os.getenv("BSKY_APP_PASSWORD")

    print(f"Logging in with handle {handle} and password {passwd}", file=sys.stderr)

    if not handle or not passwd:
        sys.stderr.write("ERROR: Please set BSKY_HANDLE and BSKY_APP_PASSWORD\n")
        sys.exit(1)

    client = Client()
    client.login(handle, passwd)
    return client


def search_hashtag_posts(client, hashtag):
    """Use the modern namespaced method: client.app.bsky.feed.search_posts."""
    cursor = None
    fetched = 0
    query = hashtag

    while True:
        resp = client.app.bsky.feed.search_posts(
            params={
                "q": query,
                "limit": PAGE_LIMIT,
                "cursor": cursor,
            }
        )

        posts = resp.posts or []
        if not posts:
            break

        for p in posts:
            if fetched >= MAX_POSTS:
                return
            text = getattr(p.record, "text", "") or ""
            if (hashtag.lower() in text.lower() or
                f"#{hashtag}".lower() in text.lower()):
                fetched += 1
                yield p

        cursor = resp.cursor
        if not cursor:
            break

        time.sleep(SLEEP)


def post_to_row(p):
    rec = p.record

    text = getattr(rec, "text", "") or ""
    created_at = getattr(rec, "created_at", "")
    uri = getattr(p, "uri", "")
    cid = getattr(p, "cid", "")

    author = getattr(p, "author", None)
    ah = getattr(author, "handle", "") if author else ""
    adisplay = getattr(author, "display_name", "") if author else ""
    adid = getattr(author, "did", "") if author else ""

    like = getattr(p, "like_count", 0) or 0
    repost = getattr(p, "repost_count", 0) or 0
    reply = getattr(p, "reply_count", 0) or 0
    quote = getattr(p, "quote_count", 0) or 0

    # Convert at:// URI → web URL
    web_url = ""
    if uri.startswith("at://"):
        try:
            _, did, _, rkey = uri.split("/", 3)
            web_url = f"https://bsky.app/profile/{did}/post/{rkey}"
        except Exception:
            pass

    return {
        "created_at": created_at,
        "author_handle": ah,
        "author_display_name": adisplay,
        "author_did": adid,
        "uri": uri,
        "cid": cid,
        "web_url": web_url,
        "text": text.replace("\n", " ").strip(),
        "like_count": like,
        "repost_count": repost,
        "reply_count": reply,
        "quote_count": quote,
    }



def main():
    client = get_client()

    print(f"Searching for #{HASHTAG} ...")
    rows = []

    for p in search_hashtag_posts(client, HASHTAG):
        rows.append(post_to_row(p))

    print(f"Found {len(rows)} matching posts")

    if not rows:
        print("No posts found.")
        return

    fieldnames = [
        "created_at", "author_handle", "author_display_name", "author_did",
        "uri", "cid", "web_url", "text",
        "like_count", "repost_count", "reply_count", "quote_count"
    ]

    with open(OUTPUT_TSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Saved to {OUTPUT_TSV}")


if __name__ == "__main__":
    main()