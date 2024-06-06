import tkinter as tk
import random
from tkinter import messagebox
from shapely.geometry import Point, Polygon
import numpy as np
from scipy.spatial import Voronoi

NUMBER = 1000


class CircleApplication(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Star Polygon and Inscribed Circle")

        self.points = []
        self.scaleFactor = 1.0
        self.lastX, self.lastY = 0, 0

        self.create_widgets()

    def create_widgets(self):
        controls = tk.Frame(self)
        controls.pack(side=tk.TOP, fill=tk.X)

        tk.Label(controls, text="n: ").pack(side=tk.LEFT)
        self.n_field = tk.Entry(controls)
        self.n_field.pack(side=tk.LEFT)

        tk.Label(controls, text="m: ").pack(side=tk.LEFT)
        self.m_field = tk.Entry(controls)
        self.m_field.pack(side=tk.LEFT)

        tk.Label(controls, text="X Coordinates: ").pack(side=tk.LEFT)
        self.x_field = tk.Entry(controls, width=20)
        self.x_field.pack(side=tk.LEFT)

        tk.Label(controls, text="Y Coordinates: ").pack(side=tk.LEFT)
        self.y_field = tk.Entry(controls, width=20)
        self.y_field.pack(side=tk.LEFT)

        self.generate_button = tk.Button(controls, text="Generate", command=self.generate)
        self.generate_button.pack(side=tk.LEFT)

        self.generate_randomly = tk.Button(controls, text="Randomly", command=self.random_points)
        self.generate_randomly.pack(side=tk.LEFT)

        self.clear_button = tk.Button(controls, text="Clear", command=self.clear)
        self.clear_button.pack(side=tk.LEFT)

        self.canvas = tk.Canvas(self, bg="white", width=800, height=600)
        self.canvas.pack(fill=tk.BOTH, expand=True)

        self.canvas.bind("<MouseWheel>", self.zoom)

    def generate(self):
        try:
            x = list(map(float, self.x_field.get().split()))
            y = list(map(float, self.y_field.get().split()))
            x_coords = np.array(x) * 100
            y_coords = np.array(y) * 100
            n = int(self.n_field.get())
            m = int(self.m_field.get())
            if len(x_coords) != len(y_coords):
                raise ValueError("X and Y coordinates must have the same length.")
        except ValueError as e:
            messagebox.showerror("Invalid input", str(e))
            return

        self.points = [Point(x, y) for x, y in zip(x_coords, y_coords)]

        if len(self.points) < n:
            messagebox.showerror("Invalid input", "Number of points must be at least n.")
            return

        self.draw_points()
        star_shape = self.form_star_shape(self.points, n, m)
        self.draw_star_shape(star_shape)
        self.draw_voronoi_diagram(star_shape)
        self.draw_largest_inscribed_circle(star_shape)

    def random_points(self):
        points = []
        for i in range(NUMBER):
            points.append(Point(random.randint(0, 100) * 10, random.randint(0, 100) * 10))
        self.points = points
        self.draw_points()
        star_shape = self.form_star_shape(self.points, NUMBER, 1)
        self.draw_star_shape(star_shape)
        self.draw_voronoi_diagram(star_shape)
        self.draw_largest_inscribed_circle(star_shape)

    def clear(self):
        self.canvas.delete("all")
        self.points = []

    def draw_points(self):
        self.canvas.delete("all")
        for point in self.points:
            self.canvas.create_oval(point.x - 2, point.y - 2, point.x + 2, point.y + 2, fill="black")

    def form_star_shape(self, points, n, m):
        if len(points) < n:
            raise ValueError("Not enough points to form the star shape.")

        center = self.find_centroid(points)

        sorted_points = sorted(points, key=lambda p: np.arctan2(p.y - center.y, p.x - center.x))

        star_shape = []
        for i in range(n):
            index = (i * m) % n
            star_shape.append(sorted_points[index])

        return star_shape

    def draw_star_shape(self, star_shape):
        for i in range(len(star_shape)):
            p1 = star_shape[i]
            p2 = star_shape[(i + 1) % len(star_shape)]
            self.canvas.create_line(p1.x, p1.y, p2.x, p2.y, fill="red", width=2)

    def draw_voronoi_diagram(self, star_shape):
        vor = Voronoi(np.array([[p.x, p.y] for p in star_shape], dtype=np.float32))
        for x, y in vor.vertices:
            self.canvas.create_oval(x - 2, y - 2, x + 2, y + 2, fill="blue")

    def find_centroid(self, points):
        centroid = Point(sum(p.x for p in points) / len(points), sum(p.y for p in points) / len(points))
        return centroid

    def draw_largest_inscribed_circle(self, star_shape):
        star_polygon = Polygon([(p.x, p.y) for p in star_shape])
        vor = Voronoi(np.array([[p.x, p.y] for p in star_shape], dtype=np.float32))

        largest_circle = None
        max_radius = 0

        for point in vor.vertices:
            point_point = Point((point[0], point[1]))
            if star_polygon.contains(point_point):
                radius = Point(point_point).distance(star_polygon.exterior)
                if radius > max_radius:
                    max_radius = radius
                    largest_circle = point_point

        if largest_circle:
            x, y = largest_circle.x, largest_circle.y
            self.canvas.create_oval(x - max_radius, y - max_radius, x + max_radius, y + max_radius, outline="green",
                                    width=2)
            messagebox.showinfo("The largest inner circle",
                                f"Coordinates: {np.round(x / 100, 2), np.round(y / 100, 2)}, Radius {np.round(max_radius / 100, 2)}")

    def zoom(self, event):
        factor = 1.1 if event.delta > 0 else 0.9
        self.scaleFactor *= factor
        self.canvas.scale("all", self.canvas.winfo_width() / 2, self.canvas.winfo_height() / 2, factor, factor)


if __name__ == "__main__":
    app = CircleApplication()
    app.mainloop()
