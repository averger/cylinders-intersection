import { defineConfig } from "vite";
import { svelte } from "@sveltejs/vite-plugin-svelte";
import tailwindcss from "@tailwindcss/vite";

const host = process.env.TAURI_DEV_HOST;

export default defineConfig({
  plugins: [svelte(), tailwindcss()],
  clearScreen: false,
  server: {
    host: host || "127.0.0.1",
    port: 5173,
    strictPort: true,
    hmr: host
      ? { protocol: "ws", host, port: 5183 }
      : undefined,
    watch: {
      ignored: ["**/src-tauri/**"],
    },
  },
  envPrefix: ["VITE_", "TAURI_ENV_*"],
  build: {
    outDir: "dist",
    emptyOutDir: true,
    target:
      process.env.TAURI_ENV_PLATFORM === "windows" ? "chrome105" : "es2022",
    sourcemap: !!process.env.TAURI_ENV_DEBUG,
    minify: !process.env.TAURI_ENV_DEBUG ? "esbuild" : false,
  },
});
