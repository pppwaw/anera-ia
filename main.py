import math
from math import sqrt
from collections import namedtuple
import api

Coll = namedtuple('Coll', ['x', 'y', "vx", "vy", "type", "mass", "time"])
cishu = 0


def cal_will_collide(x1, y1, vx1, vy1, x2, y2, vx2, vy2, r1, r2) -> (bool, float):
    if (x1 - x2) ** 2 + (y1 - y2) ** 2 == (r1 + r2) ** 2:
        return True, 0.01
    # point1_t = (x1+vx1*t, y1+vy1*t)
    # point2_t = (x2+vx2*t, y2+vy2*t)
    # point1_t and point2_t's distance is radius1+radius2
    # (x1+vx1*t-x2-vx2*t)^2+(y1+vy1*t-y2-vy2*t)^2 = (radius1+radius2)^2
    # print(x1, y1, vx1, vy1, x2, y2, vx2, vy2, r1, r2, sep=',')
    try:
        up1 = -vx1 * x1 + vx1 * x2 + vx2 * x1 - vx2 * x2 - vy1 * y1 + vy1 * y2 + vy2 * y1 - vy2 * y2
        up2 = sqrt(
            r1 ** 2 * vx1 ** 2 - 2 * r1 ** 2 * vx1 * vx2 + r1 ** 2 * vx2 ** 2 + r1 ** 2 * vy1 ** 2 - 2 * r1 ** 2 * vy1 * vy2 + r1 ** 2 * vy2 ** 2 + 2 * r1 * r2 * vx1 ** 2 - 4 * r1 * r2 * vx1 * vx2 + 2 * r1 * r2 * vx2 ** 2 + 2 * r1 * r2 * vy1 ** 2 - 4 * r1 * r2 * vy1 * vy2 + 2 * r1 * r2 * vy2 ** 2 + r2 ** 2 * vx1 ** 2 - 2 * r2 ** 2 * vx1 * vx2 + r2 ** 2 * vx2 ** 2 + r2 ** 2 * vy1 ** 2 - 2 * r2 ** 2 * vy1 * vy2 + r2 ** 2 * vy2 ** 2 - vx1 ** 2 * y1 ** 2 + 2 * vx1 ** 2 * y1 * y2 - vx1 ** 2 * y2 ** 2 + 2 * vx1 * vx2 * y1 ** 2 - 4 * vx1 * vx2 * y1 * y2 + 2 * vx1 * vx2 * y2 ** 2 + 2 * vx1 * vy1 * x1 * y1 - 2 * vx1 * vy1 * x1 * y2 - 2 * vx1 * vy1 * x2 * y1 + 2 * vx1 * vy1 * x2 * y2 - 2 * vx1 * vy2 * x1 * y1 + 2 * vx1 * vy2 * x1 * y2 + 2 * vx1 * vy2 * x2 * y1 - 2 * vx1 * vy2 * x2 * y2 - vx2 ** 2 * y1 ** 2 + 2 * vx2 ** 2 * y1 * y2 - vx2 ** 2 * y2 ** 2 - 2 * vx2 * vy1 * x1 * y1 + 2 * vx2 * vy1 * x1 * y2 + 2 * vx2 * vy1 * x2 * y1 - 2 * vx2 * vy1 * x2 * y2 + 2 * vx2 * vy2 * x1 * y1 - 2 * vx2 * vy2 * x1 * y2 - 2 * vx2 * vy2 * x2 * y1 + 2 * vx2 * vy2 * x2 * y2 - vy1 ** 2 * x1 ** 2 + 2 * vy1 ** 2 * x1 * x2 - vy1 ** 2 * x2 ** 2 + 2 * vy1 * vy2 * x1 ** 2 - 4 * vy1 * vy2 * x1 * x2 + 2 * vy1 * vy2 * x2 ** 2 - vy2 ** 2 * x1 ** 2 + 2 * vy2 ** 2 * x1 * x2 - vy2 ** 2 * x2 ** 2)
        down = (vx1 ** 2 - 2 * vx1 * vx2 + vx2 ** 2 + vy1 ** 2 - 2 * vy1 * vy2 + vy2 ** 2)
        t1 = (up1 - up2) / down

        t2 = (up1 + up2) / down

    except:
        return False, 0
    # print("cal", t1, t2)
    if math.isnan(t1) or math.isnan(t2):
        return False, 0
    if t1 < 0 and t2 < 0:
        return False, 0
    if t1 < 0:
        tt = t2
    elif t2 < 0:
        tt = t1
    else:
        tt = min(t1, t2)
    if (x1 * tt < 0 or y1 * tt < 0):
        return False, 0
    else:
        return True, tt


def cal_shouyi(mass, coll_list, pen: bool) -> float:
    if pen:

        shouyi = mass * (1 - api.SHOOT_AREA_RATIO)
    else:
        shouyi = mass
    for i in coll_list:
        if mass <= i[0].mass:
            return -mass
        shouyi += i[0].mass
    shouyi -= mass
    return shouyi


def get_coll_list(context, x, y):
    coll_list = []
    sum_col = 0
    for j in context.monsters + context.other_players + context.npc:

        cal, t = cal_will_collide(context.me.x, context.me.y, x, y, j.x, j.y, j.vx, j.vy, context.me.radius, j.radius)

        if cal:
            coll_list.append((j, t))
            sum_col += t
    coll_list.sort(key=lambda atom: atom[1])

    return coll_list, sum_col


def update(context: api.RawContext):
    global cishu
    if context.step % 10 != 0:
        return None
    print(cishu, "me", context.me.x, context.me.y, context.me.vx, context.me.vy, context.me.mass, context.me.radius)
    x, y = context.me.vx, context.me.vy
    # print(" ".join([i.type for i in context.monsters + context.other_players + context.npc]))
    coll_list, sum_col = get_coll_list(context, x, y)
    if len(coll_list) != 0:
        best_shouyi = cal_shouyi(context.me.mass, coll_list, False)
        t = sum_col / sqrt(x ** 2 + y ** 2)
        if best_shouyi < 0:
            best_shouyi -=  (100/t**2)
        else:
            best_shouyi +=  (100/t**2)
        best_shouyi = round(best_shouyi, 2)
    else:
        best_shouyi = 0
    best_rad = -114
    best_angle = 0
    print(cishu, "no coll_list",
          [Coll(i[0].x, i[0].y, i[0].vx, i[0].vy, i[0].type, i[0].mass, i[1]) for i in coll_list])
    print(cishu, "no", best_shouyi)
    angles = []
    for i in context.monsters + context.other_players + context.npc:
        angle = api.relative_angle(context.me.x, context.me.y, i.x, i.y, 0)
        angles += [(angle + 180) % 360, (angle + 90) % 360, (angle - 90) % 360]
    print(cishu, angles)
    for i in angles:
        radian = api.angle_to_radian(i)
        x, y = context.me.get_shoot_change_velocity(radian)
        print(i, "will change v", x, y)
        coll_list, sum_col = get_coll_list(context, x, y)
        if len(coll_list) != 0:
            shouyi = cal_shouyi(context.me.mass, coll_list, True)
            t = sum_col / sqrt(x ** 2 + y ** 2)
            if shouyi < 0:
                shouyi -=(100/t**2)
            else:
                shouyi +=  (100/t**2)
            shouyi = round(shouyi, 2)
        else:
            shouyi = 0
        print(cishu, i, "coll_list",
              [Coll(i[0].x, i[0].y, i[0].vx, i[0].vy, i[0].type, i[0].mass, i[1]) for i in coll_list])
        print(cishu, i, "shouyi", shouyi)
        if shouyi > best_shouyi:
            best_shouyi = shouyi
            best_rad = radian
            best_angle = i
    print(cishu, "best", best_angle, best_rad, best_shouyi)
    if best_rad != -114:
        cishu += 1
        return best_rad
