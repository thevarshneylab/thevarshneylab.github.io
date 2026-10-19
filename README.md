# thevarshneylab.github.io

Organization landing page for the Varshney Lab, served at <https://thevarshneylab.github.io/>.

The root `index.html` is a single static page listing the open-source tools from the lab. Add a new tool by copying an `<a class="tool">` block.

Tools served from subfolders of this repository:

| Path | Tool |
|---|---|
| [`pegrna-designer/`](pegrna-designer/) | [Danio pegRNA Designer](https://thevarshneylab.github.io/pegrna-designer/) — prime editing guide design for zebrafish |

Each has its own README, licence and citation file. `pegrna-designer` has a test suite run by CI on every push (`.github/workflows/test.yml`); run it locally with `node pegrna-designer/test.js`.

GitHub serves this at the organization root only while the repository is named exactly `thevarshneylab.github.io`.
