use std::net::SocketAddr;

use anyhow::Context;
use clap::Parser;
use cylix::api::router;
use tower_http::trace::TraceLayer;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt, EnvFilter};

#[derive(Parser, Debug)]
#[command(author, version, about = "Cylix — cylinder-intersection cut templates at scale 1:1")]
struct Cli {
    /// Address to bind on.
    #[arg(long, default_value = "127.0.0.1")]
    host: String,
    /// TCP port.
    #[arg(long, default_value_t = 8787)]
    port: u16,
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    tracing_subscriber::registry()
        .with(EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")))
        .with(tracing_subscriber::fmt::layer().with_target(false))
        .init();

    let args = Cli::parse();
    let addr: SocketAddr = format!("{}:{}", args.host, args.port)
        .parse()
        .with_context(|| format!("invalid host:port {}:{}", args.host, args.port))?;

    let app = router().layer(TraceLayer::new_for_http());

    tracing::info!("listening on http://{addr}");
    let listener = tokio::net::TcpListener::bind(addr).await?;
    axum::serve(listener, app).await?;
    Ok(())
}
