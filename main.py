import math
from collections import namedtuple
from math import sin, cos, sqrt
from queue import Queue
import time
import api
from api import Atom

data = namedtuple("data", ["is_space", "data"])
step = ""
start_time = 0
q = Queue()
spaces = 0
AtomTuple = namedtuple(
    "Atom", ["x", "y", "vx", "vy", "radius", "radian", "mass", "type", "id"]
)
SHANBI_CISHU = 5
TARGET_CISHU = 5
MAX_SPEED = 30
SHANBI_TIME = 2


def mirror_atoms(atoms: list[api.Atom]):
    rtn = atoms.copy()
    height = api.get_context().get_env_height()  # 地图高度
    width = api.get_context().get_env_width()  # 地图宽度
    for i in atoms:  # 四个墙壁都mirror
        rtn.append(
            AtomTuple(-i.x, i.y, -i.vx, i.vy, i.radius, api.relative_radian(0, 0, -i.vx, i.vy), i.mass, i.type, i.id))
        rtn.append(
            AtomTuple(i.x, -i.y, i.vx, -i.vy, i.radius, api.relative_radian(0, 0, i.vx, -i.vy), i.mass, i.type, i.id))
        rtn.append(
            AtomTuple(2 * width - i.x, i.y, -i.vx, i.vy, i.radius, api.relative_radian(0, 0, -i.vx, i.vy), i.mass,
                      i.type, i.id))
        rtn.append(
            AtomTuple(i.x, 2 * height - i.y, i.vx, -i.vy, i.radius, api.relative_radian(0, 0, i.vx, -i.vy), i.mass,
                      i.type, i.id))
    return rtn


def distance_to(me, i):
    return sqrt((me.x - i.x) ** 2 + (me.y - i.y) ** 2)


def qw_c(mass, t):
    return mass * 10 / t


def speeds(me: api.Atom, atom: api.Atom) -> tuple[
    float, float, float, float]:  # lianxian speed, chuizhi speed, lianxian radian, chuizhi radian
    d = api.distance(me.x, me.y, atom.x, atom.y)
    u_xiangdui = ((atom.x - me.x) / d, (atom.y - me.y) / d)  # 连线方向上的单位向量
    v_xiangdui = (me.vx - atom.vx, me.vy - atom.vy)  # 相对速度
    # 沿连线方向上的速度
    v_lianxian = v_xiangdui[0] * u_xiangdui[0] + v_xiangdui[1] * u_xiangdui[1]
    # 垂直连线方向上的速度
    u_chuizhi = (-u_xiangdui[1], u_xiangdui[0])
    v_chuizhi = v_xiangdui[0] * u_chuizhi[0] + v_xiangdui[1] * u_chuizhi[1]
    return v_lianxian, v_chuizhi, api.relative_angle(0, 0, *u_xiangdui), api.relative_angle(0, 0, *u_chuizhi)


def Angle(me: api.Atom, atom: api.Atom, cishu) -> list[float]:
    rtn = []
    d = api.distance(me.x, me.y, atom.x, atom.y)
    u_xiangdui = ((atom.x - me.x) / d, (atom.y - me.y) / d)  # 连线方向上的单位向量
    v_xiangdui = (me.vx - atom.vx, me.vy - atom.vy)  # 相对速度
    # 沿连线方向上的速度
    v_lianxian, v_chuizhi, ang_lianxian, ang_chuizhi = speeds(me, atom)
    # print(f"angle_xiangdui = {round(api.relative_angle(0, 0, *u_xiangdui),3)}, angle_v_xiangdui = {round(api.relative_angle(0, 0,* v_xiangdui),3)}")
    # print(f"v_lianxian = {v_lianxian}, v_chuizhi = {v_chuizhi}")
    ang_pen_lian = api.angle_to_radian(ang_lianxian + 180)
    ang_pen_chui = api.angle_to_radian(ang_chuizhi)
    if abs(v_chuizhi) > 0.1 and cishu > 0:
        if abs(v_chuizhi) < 10.2:
            v_shuiping = sqrt(10.2 ** 2 - v_chuizhi ** 2)
            pen_x, pen_y = v_chuizhi * cos(ang_pen_chui) + v_shuiping * cos(ang_pen_lian), v_chuizhi * sin(
                ang_pen_chui) + v_shuiping * sin(ang_pen_lian)
            rtn.append(api.relative_radian(0, 0, pen_x, pen_y))
        else:
            if v_chuizhi < 0:
                v_chuizhi = -v_chuizhi
                ang_pen_chui = api.angle_to_radian(ang_chuizhi + 180)
            cishu_chuizhi = int(v_chuizhi // 10.2)
            if cishu_chuizhi >= cishu:
                rtn = [ang_pen_chui] * cishu
            else:
                rtn = [ang_pen_chui] * cishu_chuizhi
            v_chuizhi -= cishu_chuizhi * 10.2
            if abs(v_chuizhi) > 0.1 and cishu_chuizhi < cishu:
                v_shuiping = sqrt(10.2 ** 2 - v_chuizhi ** 2)
                pen_x, pen_y = v_chuizhi * cos(ang_pen_chui) + v_shuiping * cos(ang_pen_lian), v_chuizhi * sin(
                    ang_pen_chui) + v_shuiping * sin(ang_pen_lian)
                rtn.append(api.relative_radian(0, 0, pen_x, pen_y))
    if len(rtn) < cishu:
        rtn += [ang_pen_lian] * (cishu - len(rtn))
    # print(f"angle rtn={[round(api.r2a(i),3) for i in rtn]}")
    return rtn


def shanbiAngle(me: api.Atom, other: api.Atom) -> list[float]:
    rtn = []
    dx = other.x - me.x
    dy = other.y - me.y

    # 2. 计算连线方向（弧度）
    line_direction_rad = math.atan2(dy, dx)

    # 3. 计算相对速度
    vx_relative = me.vx - other.vx
    vy_relative = me.vy - other.vy

    # 4. 计算连线方向上的速度分量
    v_along_line = vx_relative * math.cos(line_direction_rad) + vy_relative * math.sin(line_direction_rad)

    # 5. 计算需要抵消的速度向量
    vx_to_cancel = v_along_line * math.cos(line_direction_rad)
    vy_to_cancel = v_along_line * math.sin(line_direction_rad)

    # 6. 计算喷射方向（与需要抵消的速度相同）
    thrust_direction_rad = math.atan2(vy_to_cancel, vx_to_cancel)
    if len(rtn) < SHANBI_CISHU:
        rtn += [thrust_direction_rad] * (SHANBI_CISHU - len(rtn))
    print(f"angle rtn={[round(api.r2a(i), 3) for i in rtn]}")
    return rtn


def get_shoot_change_velocity(radian) -> tuple[float, float]:  # x,y
    return -cos(radian) * 10.2, -sin(radian) * 10.2


def print_atom(atom: api.Atom):
    if atom:
        return AtomTuple(
            round(atom.x, 3),
            round(atom.y, 3),
            round(atom.vx, 3),
            round(atom.vy, 3),
            round(atom.radius, 3),
            round(atom.radian, 3),
            round(atom.mass, 3),
            atom.type,
            atom.id,
        )
    else:
        return None


def hebing(angs):
    ang = []
    for i in range(len(angs[0])):
        x, y = 0, 0
        for j in angs:
            xx = cos(j[i])
            yy = sin(j[i])
            x += xx
            y += yy
        ang.append(api.relative_radian(0, 0, x, y))
    return ang


def cal_t(me: api.Atom, atom: api.Atom, vx, vy):
    vx = atom.vx - me.vx - vx
    vy = atom.vy - me.vy - vy
    dx = atom.x - me.x
    dy = atom.y - me.y
    r_sum = me.radius + atom.radius

    a = vx ** 2 + vy ** 2
    b = 2 * (dx * vx + dy * vy)
    c = dx ** 2 + dy ** 2 - r_sum ** 2

    discriminant = b ** 2 - 4 * a * c

    if discriminant < 0:
        return 999  # 无碰撞

    sqrt_discriminant = math.sqrt(discriminant)
    t1 = (-b - sqrt_discriminant) / (2 * a)
    t2 = (-b + sqrt_discriminant) / (2 * a)

    t = min(t1, t2)
    if t < 0:
        t = max(t1, t2)
    if t < 0:
        return 999  # 所有解为负，无碰撞

    return t


def whether_collide(self, other):
    # 1. 计算两个物体之间的位置差
    dx = other.x - self.x
    dy = other.y - self.y

    # 2. 计算连线方向（弧度）
    line_direction_rad = math.atan2(dy, dx)

    # 3. 计算相对速度
    vx_relative = self.vx - other.vx
    vy_relative = self.vy - other.vy

    # 4. 计算连线方向上的速度分量
    v_along_line = vx_relative * math.cos(line_direction_rad) + vy_relative * math.sin(line_direction_rad)

    # 5. 如果连线速度为负（即两物体正在接近），且距离小于两者半径之和，则认为会碰撞
    distance = math.sqrt(dx ** 2 + dy ** 2)
    if v_along_line < 0 and distance <= (self.radius + other.radius):
        return True

    return False


def handle_shanbi(context: api.RawContext):
    start_time = time.time()
    me = context.me
    enemies = [i for i in context.enemies.copy() if i.mass > me.mass]
    enemies.sort(key=lambda x: distance_to(me, x))
    # print(f"shanbi me={print_atom(me)}")
    # print(f"shanbi enemies={[print_atom(i) for i in enemies]}")
    print("******shanbi******")
    print(f"me: {print_atom(context.me)}")
    angs = []
    for e in enemies:
        if whether_collide(e, me):
            print(f"shanbi {print_atom(e)} will collide")
            t = cal_t(me, e, 0, 0)
            if t > SHANBI_TIME:
                print(f"More than {SHANBI_TIME}s, ignore")
                continue
            elif t <= -1:
                print("No collide")
                continue
            angs.append(shanbiAngle(me, e))
    print(f"angs: {angs}")

    if angs:
        ang = hebing(angs)
        print(f"final angle: {[round(api.r2a(i), 3) for i in ang]}")
        while not q.empty():
            q.get()
        for i in ang:
            q.put(data(False, i))
    print(f"handle_shanbi time: {time.time() - start_time}")
    print("******shanbi******")


def have_bigger_atom(context, me: api.Atom, i: api.Atom, cishu):
    radian = api.relative_radian(me.x, me.y, i.x, i.y)
    p_l = (
        me.x - me.radius * cos(math.pi / 2 - radian),
        me.y - me.radius * sin(math.pi / 2 - radian),
    )
    p_r = (
        me.x + me.radius * cos(math.pi / 2 - radian),
        me.y + me.radius * sin(math.pi / 2 - radian),
    )
    # print(f"have_bigger_atom radian={radian}, p_l={p_l}, p_r={p_r}")
    enemies = [
        i
        for i in context.enemies
        if i.mass >= (me.mass * (1 - api.SHOOT_AREA_RATIO) ** cishu)
    ]
    if (
            len(api.raycast(enemies, p_l, radian, distance_to(me, i)))
            + len(api.raycast(enemies, p_r, radian, distance_to(me, i)))
            > 0
    ):
        # print(f"have_bigger_atom {print_atom(i)}, continue")
        return True
    return False


def handle_target(context: api.RawContext):
    print("******target******")
    print("me: ", print_atom(context.me))
    start_time = time.time()
    max_atom, max_cishu, max_qw = None, 0, 0
    # 1. 查看不改变方向，不喷射的收益
    enemies = [i for i in context.enemies if not i.is_bullet]
    # enemies = mirror_atoms(enemies)
    direct_atoms = context.me.get_forward_direction_atoms(enemies)
    direct_atoms.sort(key=lambda x: distance_to(context.me, x))
    biggest_atom = None
    for i in direct_atoms:
        if i.mass >= context.me.mass:
            break
        if not biggest_atom or i.mass > biggest_atom.mass:
            biggest_atom = i
    if direct_atoms and biggest_atom:
        max_atom = biggest_atom
        max_cishu = 0
        max_qw = qw_c(context.me.mass + biggest_atom.mass, cal_t(context.me, biggest_atom, 0, 0))
    # 2. 查看不改变方向，喷射的收益
    if biggest_atom:
        m_cishu = int(math.log(biggest_atom.mass / context.me.mass) / math.log
        (1 - api.SHOOT_AREA_RATIO))
        x, y = 0, 0
        for i in range(m_cishu):
            me_jiaodu = api.relative_radian(0, 0, context.me.vx, context.me.vy)
            xx, yy = get_shoot_change_velocity(me_jiaodu)
            x += xx
            y += yy
            t = cal_t(context.me, biggest_atom, x, y)
            qw = qw_c(context.me.mass + biggest_atom.mass - context.me.mass * (api.SHOOT_AREA_RATIO ** m_cishu), t)
            if biggest_atom.type == "npc" or biggest_atom.type == "player":
                qw *= 0.8
            if qw > max_qw * 1.02 and t >= 0.01:
                max_qw = qw
                max_atom = biggest_atom
                max_cishu = i + 1
    # 3. 查看改变方向的收益
    enemies = [
        i for i in context.enemies if i.mass <= context.me.mass * (
                    1 - api.SHOOT_AREA_RATIO) and i.mass >= context.me.mass * api.SHOOT_AREA_RATIO and not i.is_bullet
    ]
    enemies = mirror_atoms(enemies)
    enemies.sort(key=lambda x: distance_to(context.me, x))
    # print(f"enemies: {[print_atom(i) for i in enemies]}")
    for enemy in enemies:

        cishu = int(math.log(enemy.mass / context.me.mass) / math.log(1 - api.SHOOT_AREA_RATIO))
        if cishu > TARGET_CISHU:
            cishu = TARGET_CISHU
        print(f"enemy: {print_atom(enemy)}")
        if have_bigger_atom(context, context.me, enemy, cishu):
            print(f"have_bigger_atom {print_atom(enemy)}, continue")
            continue

        angles = Angle(context.me, enemy, cishu)
        x, y = 0, 0
        for i in range(len(angles)):
            xx, yy = get_shoot_change_velocity(angles[i])
            x += xx
            y += yy
            t = cal_t(context.me, enemy, x, y)
            qw = qw_c(context.me.mass + enemy.mass - context.me.mass * (api.SHOOT_AREA_RATIO ** cishu), t)
            if enemy.type == "npc" or enemy.type == "player":
                qw *= 0.8
            if qw > max_qw * 1.02 and t >= 0.01:
                max_qw = qw
                max_atom = enemy
                max_cishu = i + 1
    if max_atom:
        print(f"final atom: {print_atom(max_atom)}")
        jd = Angle(context.me, max_atom, max_cishu)
        for i in jd:
            q.put(data(False, i))
    print(f"handle_target time: {time.time() - start_time}")
    print("******target******")


def handler(context: api.RawContext):
    # print(f"me: {print_atom(context.me)}")
    if context.step % 3 == 0:
        handle_shanbi(context)
    elif context.step % 10 == 0:
        handle_target(context)


def update(context: api.RawContext):
    # if context.step == 1:
    #     return api.a2r(75)
    if context.step % 30 == 0:
        print(f"{context.step // 30} second")
    handler(context)
    global spaces
    if spaces:
        spaces -= 1
        return None
    elif not q.empty():
        r: data = q.get()
        if r.is_space:
            spaces += r.data - 1
            return None
        else:
            return r.data
    else:
        return None
