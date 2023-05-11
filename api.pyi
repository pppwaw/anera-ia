from collections import namedtuple
from typing import List, Union, Literal, Tuple


def find_neighbors(atom: Atom, others: Union[List[Atom], None] = None, close_to_far: bool = True,
                   view_radian: float = None) -> List[Atom]: ...


def raycast(colliders: List[Atom], origin: (int, int), direction: float, max_distance: float = None) -> List[Atom]: ...


class Atom:
    being_absorbed: bool
    colliding: bool
    id: str
    is_bullet: bool
    mass: float
    position: Cartesian
    radian: Radian
    radius: float
    rotation: Radian
    type: str
    velocity: Velocity
    vx: float
    vy: float
    x: float
    y: float

    def __init__(self, id: str, type: str, x: float, y: float, vx: float, vy: float, mass: float,
                 radius: float,
                 rotation: float, being_absorbed: bool, colliding: bool): ...

    def degree_to(self, point_x: float, point_y: float) -> Degree: ...

    def degree_to_atom(self, atom: Atom) -> Degree: ...

    def distance_to(self, atom: Atom) -> float: ...

    def find_nearest_atom(self, atoms: List[Atom], ignore_bullet: bool) -> Atom: ...

    def get_atom_surface_dist(self, atom: Atom) -> float: ...

    def get_forward_direction_atoms(self, atoms: List[Atom]) -> List[Atom]: ...

    def get_shoot_change_velocity(self, theta: Union[float, Radian]) -> Tuple[float, float]: ...

    def radian_to(self, point_x: float, point_y: float) -> Radian: ...

    def radian_to_atom(self, atom: Atom) -> Radian: ...

    def shoot_atom_degree(self, atom: Atom) -> Degree: ...

    def shoot_atom_radian(self, atom: Atom) -> Radian: ...

    def shoot_point_degree(self, point_x: float, point_y: float) -> Degree: ...

    def shoot_point_radian(self, point_x: float, point_y: float) -> Radian: ...

    def whether_collide(self, atom: Atom) -> bool: ...


def a2r(angle: Union[float, Degree]) -> Radian: ...


def r2a(angle: Union[float, Radian]) -> Degree: ...


def angle_to_radian(angle: Union[float, Degree]) -> Radian: ...


def radian_to_angle(angle: Union[float, Radian]) -> Degree: ...


def relative_radian(x1: float, y1: float, x2: float, y2: float, r1: Union[float, int, None] = None) -> Radian: ...


def relative_angle(x1: float, y1: float, x2: float, y2: float, a1: Union[float, int, None] = None) -> Degree: ...


def distance(x1: float, y1: float, x2: float, y2: float) -> float: ...


SHOOT_AREA_RATIO = 0.02
GRAVITATIONAL_ACCELERATION = 9.8
DELTA_VELOCITY = 500.0


class RawContext:
    atom: List[Atom]  # 返回所有 Atom 列表
    bullets: List[Atom]  # 返回 bullet 抛射物的 Atom 列表
    done: bool  #返回本局是否结束
    enemies: List[Atom]  # 返回除context.me以外所有其它类型的Atom星体列表
    me: Atom  # 仅返回我自己的单个 Atom
    monsters: List[Atom]  # 返回 monster 类型的 Atom 列表
    npc: List[Atom]  # 返回 NPC 类型的 Atom 列表
    other_players: List[Atom]  # 返回除context.me以外的 player 类型的 Atom 列表
    sun: List[Atom]  # 返回 sun 类型的 Atom 列表
    step: int

    def get_atom_by_id(self, atom_id: str) -> Atom: ...  # 返回指定 id 对应的 Atom 星体


def get_context() -> RawContext: ...


Polar = namedtuple('Polar', ['rho', 'phi'])


class Cartesian:
    x: float
    y: float
    z: float
    t: float
    dtype: int

    def __init__(self, x: float, y: float, z: float, dtype: int): ...

    def add_polar(self, polar): ...

    def max(self) -> float: ...

    def min(self) -> float: ...

    def set(self, x: float, y: float, z: float, t: int): ...

    def sub_polar(self, polar): ...

    def to_polar(self) -> Polar: ...


class Velocity(Cartesian):
    def __init__(self, x: float, y: float = None, z: float = None, dtype: int = None): ...

    def direction(self, unit: Literal["degree", "radion"]) -> Union[Degree, Radian]: ...


class Angle(float):
    x: float

    def __init__(self, x: float): ...


class Degree(Angle):
    def add_radian(self, radian) -> Degree: ...

    def normalize(self) -> Degree: ...

    def sub_radian(self, radian) -> Degree: ...

    def to_radian(self) -> Radian: ...


class Radian(Angle):
    def add_degree(self, degree) -> Radian: ...

    def normalize(self) -> Radian: ...

    def sub_degree(self, degree) -> Radian: ...

    def to_degree(self) -> Degree: ...
