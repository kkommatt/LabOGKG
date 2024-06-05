import tkinter as tk
from tkinter import messagebox
import random
import time
from shapely.geometry import Polygon, MultiPoint
from shapely.ops import voronoi_diagram
import numpy as np


class StarShapeApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Star Shape and Voronoi Diagram")

        # Create a frame for the inputs and buttons
        self.control_frame = tk.Frame(root)
        self.control_frame.pack(side=tk.LEFT, padx=10, pady=10)

        # Create a canvas for drawing
        self.canvas = tk.Canvas(root, width=800, height=600, bg="white")
        self.canvas.pack(side=tk.RIGHT)

        # Add input fields and labels
        tk.Label(self.control_frame, text="X Coordinates:").pack()
        self.x_field = tk.Entry(self.control_frame)
        self.x_field.pack()

        tk.Label(self.control_frame, text="Y Coordinates:").pack()
        self.y_field = tk.Entry(self.control_frame)
        self.y_field.pack()

        # Add buttons
        tk.Button(self.control_frame, text="Clear", command=self.clear_fields).pack(pady=5)
        tk.Button(self.control_frame, text="Generate with Random", command=self.generate_random).pack(pady=5)
        tk.Button(self.control_frame, text="Solve Task", command=self.solve_task).pack(pady=5)

        # Add fields for N and M
        tk.Label(self.control_frame, text="N:").pack()
        self.n_field = tk.Entry(self.control_frame)
        self.n_field.pack()

        tk.Label(self.control_frame, text="M:").pack()
        self.m_field = tk.Entry(self.control_frame)
        self.m_field.pack()

        self.points = []
        self.map_statistics = {}
        self.zoom_factor = 1.5  # Zoom factor to scale the points

    def clear_fields(self):
        self.x_field.delete(0, tk.END)
        self.y_field.delete(0, tk.END)
        self.n_field.delete(0, tk.END)
        self.m_field.delete(0, tk.END)
        self.canvas.delete("all")
        self.points = []

    def generate_random(self):
        num_points = random.randint(5, 20)
        x_coords = [random.uniform(50, 750) for _ in range(num_points)]
        y_coords = [random.uniform(50, 550) for _ in range(num_points)]
        self.x_field.delete(0, tk.END)
        self.y_field.delete(0, tk.END)
        self.x_field.insert(0, " ".join(map(str, x_coords)))
        self.y_field.insert(0, " ".join(map(str, y_coords)))

    def generate_button_click(self):
        try:
            x_coords = list(map(float, self.x_field.get().split()))
            y_coords = list(map(float, self.y_field.get().split()))
            if len(x_coords) != len(y_coords):
                messagebox.showerror("Input Error", "X and Y coordinates must have the same length.")
                return
            self.points = list(zip(x_coords, y_coords))
            self.draw_points()
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers for coordinates.")

    def draw_points(self):
        self.canvas.delete("all")
        for point in self.points:
            self.draw_single_point(point)

    def draw_single_point(self, point):
        x, y = self.zoom_point(point)
        self.canvas.create_oval(x - 3, y - 3, x + 3, y + 3, fill="black")

    def zoom_point(self, point):
        x, y = point
        x = x * self.zoom_factor
        y = y * self.zoom_factor
        return x, y

    def solve_task(self):
        self.generate_button_click()
        start_time = time.time()
        try:
            n = int(self.n_field.get())
            m = int(self.m_field.get())
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid values for n and m.")
            return

        if len(self.points) < n:
            messagebox.showerror("Input Error", "Not enough points to form the star shape.")
            return

        self.canvas.delete("all")
        self.draw_points()

        star_shape = self.form_star_shape(self.points, n, m)
        self.draw_convex_hull(star_shape)

        star_polygon = Polygon(star_shape)
        voronoi_diagram_geom = voronoi_diagram(MultiPoint(self.points))

        self.draw_voronoi_diagram(voronoi_diagram_geom, star_polygon)

        end_time = time.time()
        execution_time = (end_time - start_time) * 1000  # in milliseconds
        self.map_statistics[n] = execution_time
        print(f"Execution Time: {execution_time:.2f} ms")

    def draw_voronoi_diagram(self, voronoi_diagram, star_polygon):
        for region in voronoi_diagram:
            if isinstance(region, Polygon):
                coords = list(region.exterior.coords)
                for i in range(len(coords) - 1):
                    x1, y1 = self.zoom_point(coords[i])
                    x2, y2 = self.zoom_point(coords[i + 1])
                    self.canvas.create_line(x1, y1, x2, y2, fill="orange")

    def form_star_shape(self, points, n, m):
        center = self.find_centroid(points)
        points.sort(key=lambda p: ((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2) ** 0.5)
        selected_points = points[:n]
        selected_points.sort(
            key=lambda p: (p[1] - center[1]) / (p[0] - center[0]) if p[0] != center[0] else float('inf'))

        star_shape = [selected_points[(i * m) % n] for i in range(n)]
        return star_shape

    def draw_convex_hull(self, hull):
        for i in range(len(hull)):
            p1 = self.zoom_point(hull[i])
            p2 = self.zoom_point(hull[(i + 1) % len(hull)])
            self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill="red")

    def find_centroid(self, points):
        x_coords = [p[0] for p in points]
        y_coords = [p[1] for p in points]
        centroid_x = sum(x_coords) / len(points)
        centroid_y = sum(y_coords) / len(points)
        return centroid_x, centroid_y


if __name__ == "__main__":
    root = tk.Tk()
    app = StarShapeApp(root)
    root.mainloop()
