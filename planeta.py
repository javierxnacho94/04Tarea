#!/usr/bin/env python
# -*- coding: utf-8 -*-

class Planeta(object):
    '''
    Complete el docstring.
    '''

    def __init__(self, condicion_inicial, alpha=0):
        '''
        __init__ es un método especial que se usa para inicializar las
        instancias de una clase.

        Ej. de uso:
        >> mercurio = Planeta([x0, y0, vx0, vy0])
        >> print(mercurio.alpha)
        >> 0.
        '''
        self.y_actual = condicion_inicial
        self.t_actual = 0.
        self.masa=1.
        self.alpha = alpha

    def ecuacion_de_movimiento(self):
        '''
        Implementa la ecuación de movimiento, como sistema de ecuaciónes de
        primer orden.
        '''
        x, y, vx, vy = self.y_actual
        a = self.alpha
        fx = ( (2*a / (x**2 + y**2)**2) - 1/(x**2 + y**2)**(3/2) ) * x
        fy = ( (2*a / (x**2 + y**2)**2) - 1/(x**2 + y**2)**(3/2) ) * y
        return [vx, vy, fx, fy]

    def avanza_euler(self, dt):

        x, y, vx, vy = self.y_actual
        vx, vy, fx, fy = self.ecuacion_de_movimiento()
        self.y_actual = [x + dt*vx, y + dt*vy, vx + dt*fx, vy + dt*fy]
        self.t_actual += dt

        pass

    def avanza_rk4(self, dt):

        x, y, vx, vy = self.y_actual
        t0 = self.t_actual
        vx, vy, fx, fy = self.ecuacion_de_movimiento()
        k1 = dt * [self.ecuacion_de_movimiento()[0],self.ecuacion_de_movimiento()[1],self.ecuacion_de_movimiento()[2],self.ecuacion_de_movimiento()[3]]
        self.y_actual=[x+k1[0]/2, y+k1[1]/2, vx+k1[2]/2, vy+k1[3]/2]
        self.t_actual = t0 + dt/2
        k2 = dt * [self.ecuacion_de_movimiento()[0],self.ecuacion_de_movimiento()[1],self.ecuacion_de_movimiento()[2],self.ecuacion_de_movimiento()[3]]
        self.y_actual=[x+k2[0]/2, y+k2[1]/2, vx+k2[2]/2, vy+k2[3]/2]
        self.t_actual = t0 + dt/2
        k3 = dt * [self.ecuacion_de_movimiento()[0],self.ecuacion_de_movimiento()[1],self.ecuacion_de_movimiento()[2],self.ecuacion_de_movimiento()[3]]
        self.y_actual=[x+k3[0], y+k3[1], vx+k3[2], vy+k3[3]]
        self.t_actual = t0 + dt
        k4 = dt * [self.ecuacion_de_movimiento()[0],self.ecuacion_de_movimiento()[1],self.ecuacion_de_movimiento()[2],self.ecuacion_de_movimiento()[3]]

        xi= x + (k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        yi= y + (k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        vxi= vx + (k1[2]+2*k2[2]+2*k3[2]+k4[2])/6
        vyi= vy + (k1[3]+2*k2[3]+2*k3[3]+k4[3])/6
        self.y_actual = [xi,yi,vxi,vyi]

        pass

    def avanza_verlet(self, dt):

        x, y, vx, vy = self.y_actual
        vx, vy, fx, fy = self.ecuacion_de_movimiento()

        xi= x + vx*dt + fx*(dt**2)/2
        yi= y + vy*dt + fy*(dt**2)/2
        vxi= vx + fx*dt
        vyi= vy + fy*dt
        self.y_actual = [xi,yi,vxi,vyi]
        self.t_actual += dt

        pass

    def energia_total(self):
        x, y, vx, vy = self.y_actual
        m = self.masa
        a = self.alpha

        r = (x**2 + y**2)**0.5

        E = 0.5*m*(vx**2 + vy**2) - 1/r + a/(r**2)

        return E
        
