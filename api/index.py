"""Entrada de Vercel: la API de backend/ como función serverless (vercel.json la enruta en /api/*)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "backend"))
from app import app  # noqa: E402,F401
