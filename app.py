from flask import Flask, render_template, request, make_response
import math
import re
from typing import List, Tuple

app = Flask(__name__)


@app.route("/guide", methods=["GET"])
def guide():
    return render_template("guide.html")

# =========================
# 1. Distance (Haversine)
# =========================
def haversine(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """
    Calcule la distance (en km) entre deux points (lat, lon) en degrés.
    """
    R = 6371.0  # rayon de la Terre en km
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c


# =========================
# 2. Conversion DMS → décimal
# =========================
def dms_to_decimal(deg: float, minutes: float, seconds: float, hemi: str) -> float:
    """
    Convertit DMS + hémisphère (N/S/E/W) en degrés décimaux.
    """
    value = abs(deg) + minutes / 60.0 + seconds / 3600.0
    hemi = hemi.upper()
    if hemi in ["S", "W"]:
        value = -value
    return value


def parse_single_dms(dms_str: str) -> float:
    """
    Parse une chaîne DMS ou degrés+hémisphère, par ex :
      - 12°22'52.26"N
      - 3°25'29.91"E
      - 12 22 52.26 N
      - 12.41667°N
      - 12.41667° N
    et retourne un float en degrés décimaux.
    """
    s = dms_str.strip()

    # 1) Cas simple : "12.41667° N" ou "12.41667N" etc. (degrés + hémisphère)
    simple_pattern = r"""^\s*
        (?P<deg>-?\d+(?:\.\d+)?)
        \D*
        (?P<hemi>[NSEW])
        \s*$"""
    m_simple = re.match(simple_pattern, s, re.VERBOSE | re.IGNORECASE)
    if m_simple:
        deg = float(m_simple.group("deg"))
        hemi = m_simple.group("hemi").upper()
        # minutes et secondes = 0 par défaut
        return dms_to_decimal(deg, 0.0, 0.0, hemi)

    # 2) Cas complet : 12°22'52.26"N (deg, min, sec, hemi)
    s_no_space = s.replace(" ", "")
    full_pattern = r"""^\s*
        (?P<deg>-?\d+(?:\.\d+)?)\D+
        (?P<min>\d+(?:\.\d+)?)\D+
        (?P<sec>\d+(?:\.\d+)?)\D*
        (?P<hemi>[NSEW])
        \s*$"""
    m_full = re.match(full_pattern, s_no_space, re.VERBOSE | re.IGNORECASE)
    if m_full:
        deg = float(m_full.group("deg"))
        minutes = float(m_full.group("min"))
        seconds = float(m_full.group("sec"))
        hemi = m_full.group("hemi").upper()
        return dms_to_decimal(deg, minutes, seconds, hemi)

    # 3) Cas "12 22 52.26 N" (séparé par espaces)
    parts = s.replace("°", " ").replace("'", " ").replace("\"", " ").split()
    if len(parts) == 4:
        deg = float(parts[0])
        minutes = float(parts[1])
        seconds = float(parts[2])
        hemi = parts[3].upper()
        return dms_to_decimal(deg, minutes, seconds, hemi)

    raise ValueError(f"Format DMS invalide : {dms_str}")


# =========================
# 3. Parsing latitude/longitude
# =========================
def parse_lat_lon_line(line: str) -> Tuple[float, float]:
    """
    Parse une ligne de coordonnée :
      - décimal : "12.40446, 3.38011" ou "12.40446 3.38011"
      - DMS     : "12°22'52.26\"N , 3°25'29.91\"E"
      - degrés+hémisphère : "12.41667° N, 3.51667° E"
    Retourne (lat, lon) en float (degrés décimaux).
    """
    line = line.strip()
    if not line:
        raise ValueError("Ligne vide")

    # Séparer lat/lon via virgule en priorité
    if "," in line:
        parts = [p.strip() for p in line.split(",") if p.strip()]
        if len(parts) != 2:
            raise ValueError(f"Ligne invalide (virgules) : {line}")
        lat_str, lon_str = parts
    else:
        # séparation simple par espaces (pour le décimal)
        parts = [p for p in line.split() if p]
        if len(parts) < 2:
            raise ValueError(f"Ligne invalide : {line}")
        lat_str, lon_str = parts[0], parts[1]

    def is_dms(s: str) -> bool:
        return any(ch in s for ch in ["°", "'", "\"", "N", "S", "E", "W", "n", "s", "e", "w"])

    # latitude
    if is_dms(lat_str):
        lat = parse_single_dms(lat_str)
    else:
        lat = float(lat_str)

    # longitude
    if is_dms(lon_str):
        lon = parse_single_dms(lon_str)
    else:
        lon = float(lon_str)

    return lat, lon

def dedupe_points(points, names, eps=1e-6):
    """
    Supprime les doublons de points (lat, lon) avec une tolérance eps.
    Retourne: (points_uniques, noms_uniques, doublons)
    doublons = liste de dict: {"dropped_name":..., "kept_name":..., "lat":..., "lon":...}
    """
    unique_points = []
    unique_names = []
    duplicates = []

    for (lat, lon), nm in zip(points, names):
        found = None
        for i, (ulat, ulon) in enumerate(unique_points):
            if abs(lat - ulat) <= eps and abs(lon - ulon) <= eps:
                found = i
                break

        if found is None:
            unique_points.append((lat, lon))
            unique_names.append(nm)
        else:
            duplicates.append({
                "dropped_name": nm,
                "kept_name": unique_names[found],
                "lat": lat,
                "lon": lon
            })

    return unique_points, unique_names, duplicates


def parse_points_block(block: str) -> List[Tuple[float, float]]:
    points = []
    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue
        points.append(parse_lat_lon_line(line))
    return points


# =========================
# 4. TSP heuristique (NN + 2-opt)
# =========================
def nearest_neighbor_tour(points: List[Tuple[float, float]], start_index: int = 0) -> List[int]:
    """
    Construit un tour par heuristique du plus proche voisin.
    Retourne une liste d'indices (ordre de visite, circuit fermé).
    """
    n = len(points)
    unvisited = set(range(n))
    tour = [start_index]
    unvisited.remove(start_index)

    current = start_index
    while unvisited:
        next_point = min(
            unvisited,
            key=lambda j: haversine(points[current][0], points[current][1],
                                    points[j][0], points[j][1])
        )
        tour.append(next_point)
        unvisited.remove(next_point)
        current = next_point

    # on revient au point de départ pour faire un circuit fermé
    tour.append(start_index)
    return tour


def two_opt(points: List[Tuple[float, float]], tour: List[int]) -> List[int]:
    """
    Améliore le tour avec l'algorithme 2-opt.
    """
    improved = True
    best_tour = tour[:]
    n = len(tour)

    def tour_length(t: List[int]) -> float:
        return sum(
            haversine(
                points[t[i]][0], points[t[i]][1],
                points[t[i + 1]][0], points[t[i + 1]][1]
            )
            for i in range(len(t) - 1)
        )

    best_distance = tour_length(best_tour)

    while improved:
        improved = False
        for i in range(1, n - 2):
            for k in range(i + 1, n - 1):
                new_tour = best_tour[:i] + best_tour[i:k+1][::-1] + best_tour[k+1:]
                new_distance = tour_length(new_tour)
                if new_distance < best_distance:
                    best_tour = new_tour
                    best_distance = new_distance
                    improved = True
    return best_tour




# =========================
# 5. Génération KML
# =========================
def build_kml(points: List[Tuple[float, float]],
              tour: List[int],
              names: List[str]) -> str:
    """
    Génère le contenu KML :
    - un Placemark par point
    - un LineString suivant le tour
    """
    def kml_point(name: str, lat: float, lon: float) -> str:
        return f"""
    <Placemark>
      <name>{name}</name>
      <Point>
        <coordinates>{lon},{lat},0</coordinates>
      </Point>
    </Placemark>"""

    placemarks = []
    for idx, (lat, lon) in enumerate(points):
        placemarks.append(kml_point(names[idx], lat, lon))

    coords_str = "\n        ".join(
        f"{points[i][1]},{points[i][0]},0" for i in tour
    )

    linestring = f"""
    <Placemark>
      <name>Itineraire TSP</name>
      <Style>
        <LineStyle>
          <width>3</width>
        </LineStyle>
      </Style>
      <LineString>
        <tessellate>1</tessellate>
        <coordinates>
        {coords_str}
        </coordinates>
      </LineString>
    </Placemark>"""

    kml_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>Itineraire_TSP</name>
    {''.join(placemarks)}
    {linestring}
  </Document>
</kml>
"""
    return kml_content


# =========================
# 6. Routes Flask
# =========================
@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "GET":
        return render_template("index.html")

    # POST : traitement du formulaire
    start_text = request.form.get("start_point", "").strip()
    coords_text = request.form.get("coords", "").strip()

    if not start_text or not coords_text:
        return render_template(
            "index.html",
            error="Merci de renseigner le point de départ et les autres points.",
            start_point=start_text,
            coords=coords_text
        )

    try:
        start = parse_lat_lon_line(start_text)
        other_points = parse_points_block(coords_text)
    except Exception as e:
        return render_template(
            "index.html",
            error=f"Erreur de parsing des coordonnées : {e}",
            start_point=start_text,
            coords=coords_text
        )

    points = [start] + other_points
    names = ["Start"] + [f"P{i}" for i in range(1, len(points))]
    
    points, names, duplicates = dedupe_points(points, names, eps=1e-6)

    if len(points) < 2:
        return render_template(
            "index.html",
            error="Après suppression des doublons, il ne reste pas assez de points pour calculer un itinéraire.",
            start_point=start_text,
            coords=coords_text
        )


    # TSP (NN + 2-opt)
    start_index = 0
    tour_nn = nearest_neighbor_tour(points, start_index=start_index)
    tour_opt = two_opt(points, tour_nn)

    # Construction du tableau (étapes + distances)
    route_rows = []
    total_distance = 0.0
    for step_idx, idx_point in enumerate(tour_opt):
        lat, lon = points[idx_point]
        name = names[idx_point]

        if step_idx == 0:
            leg_distance = 0.0
        else:
            prev_idx = tour_opt[step_idx - 1]
            leg_distance = haversine(
                points[prev_idx][0], points[prev_idx][1],
                lat, lon
            )
            total_distance += leg_distance

        route_rows.append({
            "step": step_idx + 1,
            "name": name,
            "lat": lat,
            "lon": lon,
            "leg_distance": leg_distance,
            "cum_distance": total_distance
        })

    return render_template(
        "index.html",
        start_point=start_text,
        coords=coords_text,
        route_rows=route_rows,
        total_distance=total_distance,
        duplicates=duplicates
    )


@app.route("/download_kml", methods=["POST"])
def download_kml():
    # On refait le calcul pour générer le KML (stateless)
    start_text = request.form.get("start_point", "").strip()
    coords_text = request.form.get("coords", "").strip()

    if not start_text or not coords_text:
        return "Coordonnées manquantes pour générer le KML", 400

    try:
        start = parse_lat_lon_line(start_text)
        other_points = parse_points_block(coords_text)
    except Exception as e:
        return f"Erreur de parsing des coordonnées : {e}", 400

    points = [start] + other_points
    names = ["Start"] + [f"P{i}" for i in range(1, len(points))]
    
    points, names, duplicates = dedupe_points(points, names, eps=1e-6)

    if len(points) < 2:
        return "Après suppression des doublons, il ne reste pas assez de points.", 400


    start_index = 0
    tour_nn = nearest_neighbor_tour(points, start_index=start_index)
    tour_opt = two_opt(points, tour_nn)

    kml_content = build_kml(points, tour_opt, names)

    response = make_response(kml_content)
    response.headers["Content-Type"] = "application/vnd.google-earth.kml+xml"
    response.headers["Content-Disposition"] = "attachment; filename=itineraire_tsp.kml"
    return response


if __name__ == "__main__":
    app.run(debug=True)
