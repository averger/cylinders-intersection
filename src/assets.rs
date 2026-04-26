//! Embeds the built Svelte SPA into the Rust binary.  The `web/dist` directory
//! is produced by `npm run build` inside `./web`.  When the directory does not
//! exist yet (e.g. during a first build of the backend in isolation), the
//! struct is still valid; the static fallback handler will then return a clear
//! message instead of serving HTML.

use rust_embed::RustEmbed;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/web/dist"]
pub struct Assets;
