# Docker Support for NetNGlyc / NetPhos Pipelines

This directory hosts the shared container image used by the NetNGlyc and NetPhos pipelines when a native binary is unavailable (Apple Silicon) or when a self-contained Linux environment is desired.

---

## 1. Contents

| File | Purpose |
|------|---------|
| `Dockerfile` | 32-bit Linux image with NetNGlyc 1.0, NetPhos, tcsh, etc. |
| `build-container.sh` | Builds `netnglyc:latest`, adding compatibility tweaks for Apple Silicon Macs. |
| `signalp_stub` | Helper invoked inside the container when host SignalP 6.0 is not available. |
| `netnglyc_stub` | Minimal stub to keep Docker workflows functional if the licensed binary is missing. |

---

## 2. Required Download

### NetNGlyc 1.0 distribution
1. Register at [DTU Bioinformatics](https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0).
2. Download `netNglyc.tar.gz` (academic use only).
3. Place the archive directly in this `docker/` folder before running `build-container.sh`.

---

## 3. Building the Image

```bash
cd BioFeatureFactory/docker
./build-container.sh
```

The script:
- Extracts `netNglyc.tar.gz` into the container.
- Installs 32-bit libraries and tcsh required by the licensed binary.
- Adds SignalP stubs so NetNGlyc can consume host predictions.

The resulting image is tagged `netnglyc:latest`.

---

## 4. When to Use Docker

| Platform                 | Docker Needed? | Notes                                                              |
|--------------------------|----------------|--------------------------------------------------------------------|
| Apple Silicon (M1/M2/M3) | **Yes**        | NetNGlyc binaries are 32-bit Intel; container provides Linux x86.  |
| Linux x86_64             | Optional       | You can run NetNGlyc natively if tcsh + 32-bit libs are installed. |
| Intel macOS              | Optional       | Use Docker if the native binary or tcsh dependencies are missing.  |

Regardless of platform, SignalP 6.0 still runs on the host; the container expects host SignalP results or uses `signalp_stub` for fallback testing.

---

## 5. Pipeline Integration

Both pipelines reference this image when `--force-docker` is set or when no native binary is detected:

- `../netnglyc/netnglyc-pipeline.py`
- `../netphos/netphos-pipeline.py`

The container provides:
- NetNGlyc 1.0 (32-bit Linux build).
- NetPhos (where licensed).
- POSIX shell utilities required by the legacy scripts.

---

## 6. Checklist

1. Download `netNglyc.tar.gz`/`netphos.tar.gz` → place in `docker/`.
2. Install Docker Desktop (macOS) or ensure Docker Engine is available (Linux).
3. Run `./build-container.sh`.
4. In pipeline commands, use the defaults (`netnglyc:latest`) or pass `--docker-image <name>` if you retagged.
5. Keep SignalP 6.0 installed on the host so the ensemble pipelines can embed SignalP context in their outputs.

---

