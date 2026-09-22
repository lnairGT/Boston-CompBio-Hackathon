"""Modal serving entrypoint for the review UI (Lane 5).

    modal deploy e2b/modal_app/web.py      # persistent URL
    modal serve  e2b/modal_app/web.py      # hot-reloading dev URL

STATUS AS SHIPPED: this file has been import-checked and its object graph builds
locally, but **it has not been deployed**. Nothing in this repository should claim a
Modal deployment until a deploy has actually run and printed a URL.

What the entrypoint needs
-------------------------
* The ``e2b`` package on the image. It is added with ``add_local_python_source`` so
  the templates and static files travel with it.
* Python dependencies: ``fastapi``, ``jinja2``, ``pydantic>=2``, ``requests``.
* **No secret is required to serve the UI itself.** The UI holds no credential, calls
  no external API from the browser, and renders no key into a template. If a spine
  adapter needs a credential it is attached here as a ``modal.Secret`` and read
  server-side only — never passed to a template or a query string.
* Concurrency: the store is per-container and in memory, so two containers do not see
  each other's runs. For a shared demo pin ``max_containers=1`` (done below) or give
  the spine a durable store and have the UI read it through a refresher.
"""

from __future__ import annotations

import os

try:
    import modal
except ModuleNotFoundError as exc:  # pragma: no cover - local dev without modal
    raise SystemExit(
        "The modal SDK is not installed in this environment. "
        "`pip install modal` to deploy; the UI itself runs fine under uvicorn without it."
    ) from exc

APP_NAME = os.environ.get("E2B_MODAL_APP_NAME", "evidence-to-binder-ui")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "fastapi[standard]==0.115.*",
        "jinja2>=3.1",
        "pydantic>=2.7",
        "requests>=2.31",
    )
    # Ships e2b/ including ui/templates and ui/static.
    .add_local_python_source("e2b")
)

app = modal.App(APP_NAME)


@app.function(
    image=image,
    # One container keeps the in-memory run store coherent for a single demo audience.
    max_containers=1,
    scaledown_window=300,
    timeout=900,
    # secrets=[modal.Secret.from_name("anthropic")],  # only if the spine agent runs here
)
@modal.concurrent(max_inputs=20)
@modal.asgi_app()
def ui():
    """Serve the review UI.

    Fixtures are preloaded so every interface state is walkable immediately. Set
    ``E2B_UI_PRELOAD_FIXTURES=0`` to serve a clean instance with no fixture runs.
    """
    from e2b.ui.app import create_app

    return create_app()


@app.local_entrypoint()
def main() -> None:
    print(
        f"App '{APP_NAME}' defined. Deploy with:\n"
        f"    modal deploy e2b/modal_app/web.py\n"
        f"Modal prints the URL on success; do not quote a URL that was never printed."
    )
