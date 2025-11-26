from manim import *
import numpy as np


# --------------------------------------------------
# Definições de funções auxiliares para as animações
# --------------------------------------------------


# Plano tangente
def tangent_plane(point):
    """Retorna o plano tangente às subvariedades integrais que passa pelo ponto point.
    Tal ponto deve ser um np.array"""

    x0, y0, z0 = point
    v1 = np.array([1, 0, y0])  # ∂/∂x
    v2 = np.array([0, 1, x0])  # ∂/∂y

    return Surface(
        lambda u0, v0: point + u0 * v1 + v0 * v2,
        u_range=[-5, 5],
        v_range=[-5, 5],
        resolution=(6, 6),
        fill_color=GREEN,
        fill_opacity=0.7,
    )


# Updater do plano
def update_plane(plano):
    """Updater necessário para que o plano passe sempre pelo ponto dot.
    Obs: ponto dot é o mesmo ponto utilizado em todas as animações"""

    new_point = dot.get_center()
    plano.become(tangent_plane(new_point))


def parametrizacao_01(t):
    """Vai do ponto (1, -1, 1) até o ponto (-1, 1, 1) sobre a curva de nível z = x * y + 2."""

    x0 = 1 - 2 * t
    y0 = -1 + 2 * t
    z0 = x0 * y0 + 2
    return np.array([x0, y0, z0])


def parametrizacao_02(t):
    """Vai do ponto (-1, 1, 1) até o ponto (1, 1, 3) sobre a curva de nível z = x * y + 2."""

    x0 = -1 + 2 * t
    y0 = 1
    z0 = x0 * y0 + 2
    return np.array([x0, y0, z0])


def parametrizacao_03(t):
    """Vai do ponto (1, 1, 3) até o ponto (-1, -1, 3) sobre a curva de nível z = x * y + 2."""

    x0 = 1 - 2 * t
    y0 = 1 - 2 * t
    z0 = x0 * y0 + 2
    return np.array([x0, y0, z0])


def parametrizacao_04(t):
    """Vai do ponto (-1, -1, 3) até o ponto (1, 2, 0) em linha reta."""

    x0 = -1 + 2 * t
    y0 = -1 + 3 * t
    z0 = 3 - 3 * t
    return np.array([x0, y0, z0])


def subvariedade_integral(point, x_range=(-6, 6), y_range=(-6, 6), res=30):
    """Retorna a superfície de nível de F(x,y,z) = z - x * y que passa pelo ponto fornecido.
    point: array/list/np.array de formato (3,) representando (x0, y0, z0)"""

    x0, y0, z0 = point
    c = z0 - x0 * y0  # constante que define a curva de nível

    return Surface(
        lambda u0, v0: np.array([
            u0,
            v0,
            u0 * v0 + c  # z = x * y + c
        ]),
        u_range=[x_range[0], x_range[1]],
        v_range=[y_range[0], y_range[1]],
        resolution=(res, res),
        fill_opacity=1.0,
        stroke_color=WHITE,
    )


# Updater da subvariedade integral
def update_subvariedade_integral(superficie):
    """Updater necessário para que a superfície passe sempre pelo ponto dot.
    Obs: ponto dot é o mesmo ponto utilizado em todas as animações"""

    new_point = dot.get_center()
    superficie.become(subvariedade_integral(new_point))


def v(point):
    """Retorna o vetor do campo V no ponto point (como Arrow3D)."""
    """Ponto precisa ser um Dot3D"""

    x0, y0, z0 = point.get_center()
    v_p = Arrow3D(
        start=point.get_center(),
        end=point.get_center() + [1, 0, y0],
        color=ORANGE,
        resolution=8,
        stroke_width=1,
    )

    return v_p


def w(point):
    """Retorna o vetor do campo W no ponto point (como Arrow3D)."""
    """Ponto precisa ser um Dot3D"""

    x0, y0, z0 = point.get_center()
    w_p = Arrow3D(
        start=point.get_center(),
        end=point.get_center() + [0, 1, x0],
        color=PURPLE,
        resolution=8,
        stroke_width=1,
    )

    return w_p


def v_proj(point):
    """Projeta o vetor do campo V no plano xy"""

    x0, y0, z0 = point.get_center()
    return Arrow3D(
        start=np.array([x0, y0, 0]),
        end=np.array([x0 + 1, y0, 0]),
        color=ORANGE,
        resolution=8,
        stroke_width=1
    )


def w_proj(point):
    """Projeta o vetor do campo W no plano xy"""

    x0, y0, z0 = point.get_center()
    return Arrow3D(
        start=np.array([x0, y0, 0]),
        end=np.array([x0, y0 + 1, 0]),
        color=PURPLE,
        resolution=8,
        stroke_width=1
    )


def fluxo_v(t, x0, y0, z0):
    """Dado um certo t e um ponto p = (x0, y0, z0), retorna o ponto phi_t(p) (fluxo de V)."""

    return np.array([
        x0 + t,
        y0,
        z0 + t * y0
    ])


def fluxo_w(t, x0, y0, z0):
    """Dado um certo t e um ponto p = (x0, y0, z0), retorna o ponto phi_t(p) (fluxo de W)."""

    return np.array([
        x0,
        y0 + t,
        z0 + t * x0
    ])


def big_phi(u0, v0, w0):
    """Sistema de coordenadas final do teorema de Frobenius. Retorna o ponto correspondente a:
    Andar z unidades no eixo z, andar y unidades no fluxo de W, andar x unidades no fluxo de V.
    Phi (u, v, w) = phi^V_u . phi^W_v (0, 0, w)"""

    return np.array([u0, v0, w0 + u0 * v0])  # Específico para o exemplo, mudar para generalizar


def update_cod(ponto):
    """Updater necessário para que dot_cod seja sempre atualizado para ser a imagem de dot_dom pela big_phi"""

    u0, v0, w0 = axes_dom.p2c(dot_dom.get_center())  # pegar coordenadas em eixos do domínio
    x0, y0, z0 = big_phi(u0, v0, w0)  # Aplicar função big_phi
    ponto.move_to(axes_cod.c2p(x0, y0, z0))  # Mover o ponto (que será dot_cod) para a imagem de dot_dom


# ----------------------------------------------
# Criação dos elementos utilizados nas animações
# ----------------------------------------------


# Axes (eixos) principais a serem utilizados em todas as animações, exceto SistemaDeCoordenadas
axes = ThreeDAxes(
    x_range=[-6, 6, 1],
    y_range=[-6, 6, 1],
    z_range=[-6, 6, 1],
    x_length=12,
    y_length=12,
    z_length=12,
)

# Axes a serem utilizados na animação SistemaDeCoordenadas
axes_dom = ThreeDAxes(  # axes_dom é o domínio da função big_phi
    x_range=[-5, 5, 1],
    y_range=[-5, 5, 1],
    z_range=[-5, 5, 1],
    x_length=9,
    y_length=9,
    z_length=9,
).shift(RIGHT * 2.5 + DOWN * 2.5)

axes_cod = ThreeDAxes(  # axes_cod é o contradomínio da função big_phi
    x_range=[-5, 5, 1],
    y_range=[-5, 5, 1],
    z_range=[-5, 5, 1],
    x_length=9,
    y_length=9,
    z_length=9,
).shift(LEFT * 2.5 + UP * 2.5)

# Ponto que será utilizado nas animações, exceto SistemaDeCoordenadas
dot = Dot3D(np.array([-1, -1, -1]), radius=0.08, color=YELLOW)

# Pontos utilizados na animação SistemaDeCoordenadas
dot_dom = Dot3D(axes_dom.c2p(0, 0, 0), color=YELLOW)  # dom = domínio da big_phi
dot_cod = Dot3D(axes_cod.c2p(0, 0, 0), color=YELLOW)  # cod = contradomínio da big_phi
dot_cod.add_updater(update_cod)                       # Updater no dot_cod

# Campos de vetores X e Y iniciais do exemplo
xs = np.linspace(-3, 3, 4)  # Coordenadas dos pontos que serão usados. No total, são 64 vetores (4 ** 3)
ys = np.linspace(-3, 3, 4)
zs = np.linspace(-3, 3, 4)

vector_field_X = []  # Cria as listas que serão utilizadas para guardar os campos de vetores
vector_field_Y = []

# Geração completa dos campos X e Y
for x in xs:
    for y in ys:
        for z in zs:
            p = np.array([x, y, z])

            Xv = np.array([x, 1, x * (y + 1)])
            Yv = np.array([1, 0, y])
            norm_Xv = np.sqrt(x ** 2 + 1 + (x * (y + 1)) ** 2)
            norm_Yv = np.sqrt(1 + y ** 2)

            arrowX = Arrow3D(
                start=p,
                end=p + 0.5 * Xv / norm_Xv,
                color=RED,
                resolution=8,
                stroke_width=1,
            )
            arrowY = Arrow3D(
                start=p,
                end=p + 0.5 * Yv / norm_Yv,
                color=GREEN,
                resolution=8,
                stroke_width=1,
            )

            vector_field_X.append(arrowX)
            vector_field_Y.append(arrowY)

# "Mini" campo de vetores X. Dos 64 vetores de vector_field_X, pega apenas 5 para mostrar como a distribuição
# é gerada por X e Y. Mesma ideia para mini_vector_field_Y
mini_vector_field_X = [vector_field_X[21],
                       vector_field_X[25],
                       vector_field_X[26],
                       vector_field_X[41],
                       vector_field_X[38]]

mini_vector_field_Y = [vector_field_Y[21],
                       vector_field_Y[25],
                       vector_field_Y[26],
                       vector_field_Y[41],
                       vector_field_Y[38]]

# Caminhos a serem utilizados para que o ponto dot caminhe sobre as subvariedades integrais
level_path_01 = ParametricFunction(
    lambda t: parametrizacao_01(t),
    t_range=(0, 1),
    color=YELLOW,
    stroke_width=2,
)
level_path_02 = ParametricFunction(
    lambda t: parametrizacao_02(t),
    t_range=(0, 1),
    color=YELLOW,
    stroke_width=2,
)
level_path_03 = ParametricFunction(
    lambda t: parametrizacao_03(t),
    t_range=(0, 1),
    color=YELLOW,
    stroke_width=2,
)
path_04 = ParametricFunction(
    lambda t: parametrizacao_04(t),
    t_range=(0, 1),
    color=YELLOW,
    stroke_width=2,
)

# Coordenadas pelas quais o ponto passará na primeira animação
start = np.array([-1, -1, -1])
stop1 = np.array([-1, 1, -1])
stop2 = np.array([-1, 1, 1])
stop3 = np.array([1, 1, -1])
end = np.array([1, -1, 1])
lista_de_pontos = [start, stop1, stop2, stop3, end]

# Criação do plano
plane = tangent_plane(start)
plane.add_updater(update_plane)  # Updater para seguir o ponto dot

# Criação das subvariedades integrais
curvas_nivel_subv_integraveis = subvariedade_integral(dot.get_center())
curvas_nivel_subv_integraveis.add_updater(update_subvariedade_integral)  # Updater para seguir o ponto dot

# Elementos para a animação SistemaDeCoordenadas:

# Linha azul no domínio, representando coordenada z do dot_dom
linha_z_dom = always_redraw(
    lambda: Line3D(
        start=axes_dom.c2p(
            0,
            0,
            0
        ),
        end=axes_dom.c2p(
            0,
            0,
            axes_dom.p2c(dot_dom.get_center())[2]
        ),
        color=BLUE,
        thickness=0.03,
    )
)

# Linha verde no domínio, representando coordenada y do dot_dom
linha_y_dom = always_redraw(
    lambda: Line3D(
        start=axes_dom.c2p(
            0,
            0,
            axes_dom.p2c(dot_dom.get_center())[2]
        ),
        end=axes_dom.c2p(
            0,
            axes_dom.p2c(dot_dom.get_center())[1],
            axes_dom.p2c(dot_dom.get_center())[2]
        ),
        color=GREEN,
        thickness=0.03,
    )
)

# Linha vermelha no domínio, representando coordenada x do dot_dom
linha_x_dom = always_redraw(
    lambda: Line3D(
        start=axes_dom.c2p(
            0,
            axes_dom.p2c(dot_dom.get_center())[1],
            axes_dom.p2c(dot_dom.get_center())[2]
        ),
        end=axes_dom.c2p(
            axes_dom.p2c(dot_dom.get_center())[0],
            axes_dom.p2c(dot_dom.get_center())[1],
            axes_dom.p2c(dot_dom.get_center())[2]
        ),
        color=RED,
        thickness=0.03,
    )
)

# OBS: primeiro é traçada a linha AZUL, depois a VERDE e por último a VERMELHA.
# Essa ordem é a que mais se adequa à animação, tendo em vista que
# iremos manter a coordenada z dos pontos constantes.

# Linha azul no contradomínio, representando coordenada z do dot_cod
linha_z_cod = always_redraw(
    lambda: Line3D(
        start=axes_cod.c2p(
            0,
            0,
            0
        ),
        end=axes_cod.c2p(
            0,
            0,
            axes_cod.p2c(dot_dom.get_center())[2]
        ),
        color=BLUE,
        thickness=0.03,
    )
)

# Linha verde no contradomínio, representando caminhar sobre o fluxo de W
linha_y_cod = always_redraw(
    lambda: Line3D(
        start=axes_cod.c2p(
            0,
            0,
            axes_cod.p2c(dot_dom.get_center())[2]
        ),
        end=axes_cod.c2p(
            *fluxo_w(
                axes_dom.p2c(dot_dom.get_center())[1],  # parâmetro t
                0,
                0,
                axes_cod.p2c(dot_dom.get_center())[2]
            )
        ),
        color=GREEN,
        thickness=0.03,
    )
)

# Linha vermelha no contradomínio, representando caminhar sobre o fluxo de V
linha_x_cod = always_redraw(
    lambda: Line3D(
        start=axes_cod.c2p(  # Note que end=(0, v, w) = fluxo_w(v, 0, 0, w)
            *fluxo_w(
                axes_dom.p2c(dot_dom.get_center())[1],  # parâmetro t
                0,
                0,
                axes_cod.p2c(dot_dom.get_center())[2]
            )
        ),
        end=axes_cod.c2p(
            *fluxo_v(
                axes_dom.p2c(dot_dom.get_center())[0],  # parâmetro t
                *fluxo_w(
                    axes_dom.p2c(dot_dom.get_center())[1],  # parâmetro t
                    0,
                    0,
                    axes_cod.p2c(dot_dom.get_center())[2]
                )
            )
        ),
        color=RED,
        thickness=0.03,
    )
)

lista_de_pontos_sistemadecoordenadas = [
    (2, -2, 0),
    (1.5, -2, 0),
    (1.5, 2, 0),
    (1, 2, 0),
    (1, -2, 0),
    (0.5, -2, 0),
    (0.5, 2, 0),
    (0, 2, 0),
    (0, -2, 0),
    (-0.5, -2, 0),
    (-0.5, 2, 0),
    (-1, 2, 0),
    (-1, -2, 0),
    (-1.5, -2, 0),
    (-1.5, 2, 0),
    (-2, 2, 0),
    (-2, -2, 0),
]


# ---------
# Animações
# ---------


class DistribuicaoCena01(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Ponto vagando pelo R^3 mostrando a distribuição em diferentes pontos
        """

        self.add(axes)
        self.add(dot, plane)

        # Movimento em linhas retas do ponto, mostrando a distribuição
        for ponto in lista_de_pontos:
            self.play(dot.animate.move_to(ponto), run_time=2, rate_func=smooth)
            self.wait()

        self.wait(5)
        self.remove(dot, plane, axes)
        self.wait()


class DistribuicaoCena02(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Ilustração dos dois campos de vetores X e Y que geram a distribuição
        """

        # Texto do campo X
        texto_campo_x = MathTex(
            r"X = x\,\frac{\partial}{\partial x} + \frac{\partial}{\partial y} + x(y+1)\,\frac{\partial}{\partial z}",
            color=RED
        ).scale(0.7).to_corner(UL)

        # Texto do campo Y
        texto_campo_y = MathTex(
            r"Y = \frac{\partial}{\partial x} + y\,\frac{\partial}{\partial z}",
            color=GREEN
        ).scale(0.7).to_corner(UL)

        self.add(axes)

        # Mostra o campo de vetores X
        axes.add(*vector_field_X)
        self.wait(2)

        self.add_fixed_in_frame_mobjects(texto_campo_x)
        self.play(Write(texto_campo_x))
        self.wait(2)

        # Rotaciona o axes para vermos melhor o campo X
        self.play(
            Rotate(
                axes,
                angle=2 * PI,
                axis=OUT,
                run_time=25,
                rate_func=linear
            )
        )

        self.wait(2)

        axes.remove(*vector_field_X)
        self.remove(texto_campo_x)
        self.wait(2)

        # Mostra o campo de vetores Y
        self.add(*vector_field_Y)
        self.wait(2)

        self.add_fixed_in_frame_mobjects(texto_campo_y)
        self.play(Write(texto_campo_y))
        self.wait(2)

        self.move_camera(
            phi=90 * DEGREES,
            theta=0 * DEGREES,
            run_time=2,
            rate_func=smooth
        )
        self.move_camera(
            phi=90 * DEGREES,
            theta=360 * DEGREES,
            run_time=25,
            rate_func=linear
        )
        self.move_camera(
            phi=60 * DEGREES,
            theta=405 * DEGREES,
            run_time=2,
            rate_func=smooth
        )

        self.wait(2)
        self.remove(*vector_field_Y, texto_campo_y)
        self.wait(2)


class DistribuicaoCena03(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Mostra como os campos geram a distribuição
        """

        self.add(axes)

        # Mostra os campos por onde o ponto passa
        self.add(*mini_vector_field_X, *mini_vector_field_Y)

        self.add(dot, plane)
        self.wait()

        # Movimentos em linha reta
        for ponto in lista_de_pontos:
            self.play(dot.animate.move_to(ponto), run_time=2, rate_func=smooth)
            self.wait(3)
        self.wait(5)

        self.remove(axes, dot, plane, *mini_vector_field_X, *mini_vector_field_Y)
        self.wait()


class LemaDistribuicaoInvolutiva(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Texto na tela, demonstração do lema para distribuição involutiva
        """

        enunciado = Tex(
            r"\begin{flushleft}"
            r"Lema: Sejam $X, Y \in \mathfrak{X}(M)$ e "
            r"$D = \operatorname{ger}\{X,Y\}$ a distribuição gerada por tais campos. "
            r"Se $[X,Y](p) \in D(p)$ para todo ponto $p \in M$, "
            r"então a distribuição $D$ é involutiva."
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).to_corner(UP)

        dem1 = Tex(
            r"\begin{flushleft}"
            r"Demonstração: sejam $V, W \in \mathfrak{X}(M)$ campos subordinados a $D$, isto é, "
            r"$V(p), W(p) \in D(p)$ para todo $p \in M$. Mostremos que $\left[ V, W \right](p) \in D(p)$, "
            r"para todo ponto $p \in M$."
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).move_to(UP * 1.25)

        dem2 = Tex(
            r"\begin{flushleft}"
            r"Sabemos que $V = f_1 X + g_1 Y$ e $W = f_2 X + g_2 Y$, para "
            r"certas $f_1, f_2, g_1, g_2 \in C^{\infty} (M)$. Desta forma, usando a bilinearidade "
            r"do colchete de Lie, temos que "
            r"\end{flushleft}"
        ).scale(0.7).move_to(DOWN * 0.5)

        dem3 = MathTex(
            r"\begin{aligned}"
            r"[V, W]"
            r"&= [f_1 X + g_1 Y, f_2 X + g_2 Y] \\"
            r"&= [f_1 X, f_2 X] + [f_1 X, g_2 Y] + [g_1 Y, f_2 X] + [g_1 Y, g_2 Y]"
            r"\end{aligned}",
            color=WHITE
        ).scale(0.7).move_to(DOWN * 2.75)

        dem4 = Tex(
            r"Pela identidade $[fX, gY] = fg[X,Y] + f(Xg)Y - g(Yf)X$, obtemos:",
            color=WHITE
        ).scale(0.7).to_corner(UP)

        dem5 = MathTex(
            r"[f_1 X, f_2 X] = f_1 f_2 [X,X] + f_1(X f_2)X - f_2(X f_1)X, \\[6pt]"
            r"[f_1 X, g_2 Y] = f_1 g_2 [X,Y] + f_1(X g_2)Y - g_2(Y f_1)X, \\[6pt]"
            r"[g_1 Y, f_2 X] = g_1 f_2 [Y,X] + g_1(Y f_2)X - f_2(X g_1)Y, \\[6pt]"
            r"[g_1 Y, g_2 Y] = g_1 g_2 [Y,Y] + g_1(Y g_2)Y - g_2(Y g_1)Y.",
            color=WHITE
        ).scale(0.7).move_to(UP * 1.25)

        dem6 = Tex(
            r"Assim, temos $[V, W] = \alpha X + \beta Y + \theta [X, Y]$, com ",
            color=WHITE
        ).scale(0.7).move_to(DOWN)

        dem7 = Tex(
            r"$\alpha = f_1 X(f_2) - f_2 X(f_1) - g_2 Y(f_1) + g_1 Y(f_2)$;\\[6pt]"
            r"$\beta  = f_1 X(g_2) - f_2 X(g_1) + g_1 Y(f_2) - g_2 Y(g_1)$;\\[6pt]"
            r"$\theta = f_1 g_2 - g_1 f_2$.",
            color=WHITE
        ).scale(0.7).move_to(DOWN * 2.75)

        dem8 = Tex(
            r"Logo, para todo $p \in M$, \\"
            r"$[V, W](p) = \alpha(p)X(p) + \beta(p)Y(p) + \theta(p)[X, Y](p) \in D(p)$. \\"
            r"Portanto, D é involutiva.",
            color=WHITE
        ).scale(0.7).to_corner(UP)

        textos = [enunciado, dem1, dem2, dem3, dem4, dem5, dem6, dem7, dem8]

        textos_na_tela = []  # Guarda quais textos estão sendo exibidos
        for texto in textos:

            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto, run_time=3))
            self.wait(6)

            textos_na_tela.append(texto)

            if len(textos_na_tela) == 4:  # Máximo de textos por vez
                self.remove(*textos_na_tela)
                textos_na_tela.clear()

        self.remove(*textos_na_tela)
        self.wait()


class CalculoColchete(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Texto na tela, cálculo do colchete de X e Y
        """

        texto_1 = Tex(
            r"\begin{flushleft}"
            r"Utilizando o lema anterior, vamos verificar que nossa distribuição $D$ é involutiva. "
            r"Para isso, calculamos o colchete de Lie entre X e Y:"
            r"\end{flushleft}"
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"\begin{aligned}"
            r"[X, Y] "
            r"&= \sum _{j=1} ^3 (X(b_j) - Y(a_j)) \frac{\partial}{\partial u_j}\\"
            r"&= (X(1) - Y(x)) \frac{\partial}{\partial x} +"
            r"(X(0) - Y(1)) \frac{\partial}{\partial y} + "
            r"(X(y) - Y(x(y + 1))) \frac{\partial}{\partial z}\\"
            r"&= (0 - 1) \frac{\partial}{\partial x} +"
            r"(0 - 0) \frac{\partial}{\partial y} + "
            r"(1 - (y + 1)) \frac{\partial}{\partial z}\\"
            r"&= - \frac{\partial}{\partial x} +"
            r" - y \frac{\partial}{\partial z}\\"
            r"&= - Y"
            r"\end{aligned}"
        ).scale(0.7).move_to(UP * 0)

        texto_3 = Tex(
            r"\begin{flushleft}"
            r"Deste modo, sendo $p \in M$ qualquer, $[X, Y](p) = -Y(p) \in D$ e, pelo lema, $D$ é involutiva."
            r"\end{flushleft}"
        ).scale(0.7).move_to(DOWN * 3)

        textos = [texto_1, texto_2, texto_3]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto))
            self.wait(4)
        self.wait(2)

        self.remove(*textos)
        self.wait()


class TextoTeoremaFrobeniusGeometrico(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Texto na tela, enunciação do teorema geométrico de Frobenius
        """

        texto_1 = Tex(
            r"Teorema de Frobenius: uma distribuição $D$ de posto $k$ é integrável se, e somente se, for involutiva."
        ).scale(0.8).to_corner(UP)

        texto_2 = Tex(
            r"\begin{flushleft}"
            r"Para cada $p \in M$, existe um sistema cúbico de coordenadas $\varphi: U \longrightarrow \varphi(U) ="
            r"(- \varepsilon , \varepsilon)^m$, com $\varphi(p) = 0$ tal que, para todo "
            r"$a = (a_1, \dots, a_{m-k}) \in (- \varepsilon , \varepsilon)^{m-k}$, as subvariedades "
            r"$S_a = (\pi \circ \varphi)^{-1} (a)$ são subvariedades integrais de $D$ em $U$, em que "
            r"$\pi : \mathbb{R}^m \longrightarrow \mathbb{R}^{m-k}$ denota a projeção canônica. "
            r"Além disso, toda subvariedade integral conexa de $D$ em $U$ está contida em alguma $S_a$."
            r"\end{flushleft}"
        ).scale(0.7).move_to(UP)

        texto_3 = Tex(
            r"\begin{flushleft}"
            r"O teorema possui uma demonstração construtiva e seu próprio enunciado remete a essa construção"
            r"\end{flushleft}"
        ).scale(0.7).move_to(DOWN * 2)

        textos = [texto_1, texto_2, texto_3]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto), run_time=6)
            self.wait(4)
        self.wait(2)

        self.remove(*textos)
        self.wait()


class SubvariedadesIntegrais(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)
        dot.move_to(end)

        """
        CENA: Mostrar subvariedades integrais como curvas de nível
        """

        self.wait()
        self.add(axes, dot, plane)
        self.wait(2)
        self.move_camera(
            phi=60 * DEGREES,
            theta=120 * DEGREES,
            run_time=1,
            rate_func=smooth
        )
        self.wait(2)

        self.add(curvas_nivel_subv_integraveis)

        # Ponto dot percorrendo trajetórias e ilustrando distribuição tangente às subvariedades
        for curve in [level_path_01, level_path_02, level_path_03, path_04]:
            self.play(MoveAlongPath(dot, curve), run_time=4, rate_func=smooth)
            self.wait(2)
        self.wait(4)

        foliacoes = []  # Guarda 4 superfícies de nível para mostrar as foliações da distribuição
        for z0 in [-2, 0, 2, 4]:
            point = np.array([1, 2, z0])
            foliacoes.append(
                subvariedade_integral(point, x_range=(-1, 1), y_range=(-1, 1), res=10)
            )

        self.remove(plane, curvas_nivel_subv_integraveis)

        for folha in foliacoes:
            self.play(Create(folha))
            self.wait(1)
        self.wait(3)

        # Movimento de câmera para mostrar curvas de nível com mais clareza
        self.move_camera(
            phi=90 * DEGREES,
            theta=0 * DEGREES,
            run_time=3,
            rate_func=smooth
        )
        self.move_camera(
            phi=90 * DEGREES,
            theta=360 * DEGREES,
            run_time=25,
            rate_func=linear
        )
        self.move_camera(
            phi=60 * DEGREES,
            theta=45 * DEGREES,
            run_time=3,
            rate_func=smooth
        )

        self.wait()
        self.remove(dot, curvas_nivel_subv_integraveis, *foliacoes, axes)
        self.wait(5)


class TextoLema01TeoremaFrobenius(ThreeDScene):
    def construct(self):

        """
        CENA: Texto na tela, enunciação do Lema 01 teorema de Frobenius geométrico
        """

        frobenius_lema1_1 = Tex(
            r"\begin{flushleft}"
            r"Lema 1: Seja $D$ uma distribuição involutiva de posto $k$ em $M$. \\"
            r"Então, para cada $p \in M$, existem um aberto $V$ contendo $p$ e "
            r"$X_1, \dots, X_k \in \mathfrak{X}(M)$ tais que "
            r"$X_1(q), \dots, X_k(q)$ geram $D(q)$ para todo $q \in V$ e "
            r"$\left[ X_i, X_j \right] = 0$ para quaisquer $1 \leq i,j \leq k$."
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).to_corner(UP)

        frobenius_lema1_2 = Tex(
            r"\begin{flushleft}"
            r"Isto é, toda distribuição involutiva de posto $k$ pode ser localmente gerada "
            r"por $k$ campos que comutam."
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).move_to(UP)

        frobenius_lema1_3 = Tex(
            r"\begin{flushleft}"
            r"Nosso objetivo: encontrar quais são esses dois campos"
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).move_to(UP * 0)

        cena_lema_frobenius_1 = [frobenius_lema1_1, frobenius_lema1_2, frobenius_lema1_3]

        for texto in cena_lema_frobenius_1:
            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto, run_time=4))
            self.wait(2)
        self.wait(2)

        self.remove(*cena_lema_frobenius_1)
        self.wait()


class CamposVeW(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Campos V e W pi-relacionados com d/dx e d/dy
        """

        self.add(axes, dot, plane)
        self.play(dot.animate.move_to(np.array([1, 1, 1])), run_time=1, rate_func=smooth)
        self.wait()
        self.move_camera(
            phi=60 * DEGREES,
            theta=225 * DEGREES,
            run_time=1,
            rate_func=smooth
        )
        self.wait()

        # Criação dos vetores do campo V e W no ponto (1, 1, 1)
        v_original = v(dot)
        w_original = w(dot)
        self.add(v_original, w_original)
        self.add(v(dot), w(dot))
        self.wait()

        v_xy = v_proj(dot)
        w_xy = w_proj(dot)

        # Animação dos vetores v_p e w_p sendo projetados no plano xy
        self.play(
            Transform(v_original, v_xy),
            Transform(w_original, w_xy),
            run_time=3,
            rate_func=linear
        )
        self.wait()

        # Movimento de câmera para mais clareza
        self.move_camera(
            phi=60 * DEGREES,
            theta=45 * DEGREES,
            run_time=8,
            rate_func=linear
        )
        self.wait()

        # Mostrar mais vetores que também são projetados onde v_p é projetado
        pontos_ponta_vetores = [
            np.array([2, 1, -1]),
            np.array([2, 1, 0]),
            np.array([2, 1, 1]),
            np.array([2, 1, 3]),
        ]

        vetores = []
        for ponto in pontos_ponta_vetores:
            vetor = Arrow3D(
                start=np.array([1, 1, 1]),
                end=ponto,
                color=YELLOW,
                resolution=8,
                stroke_width=1,
            )

            vetores.append(vetor)

        for vetor in vetores:
            self.add(vetor)
            self.wait()
        self.wait(2)

        # Mostra a reta perpendicular ao plano xy projetada no ponto (2, 1, 0)
        reta = Line3D(
            start=np.array([2, 1, -6]),
            end=np.array([2, 1, 6]),
            color=YELLOW,
            thickness=0.01
        )
        self.play(Create(reta))
        self.wait(6)

        # Equação "diferencial x V_p = e_1" para aparecer na tela
        texto_campo_v = Tex(
            r"$d \pi (p) \cdot V_p = \frac{\partial}{\partial x} \right| _p = e_1$, para todo $p \in \mathbb{R}^3$",
            color=WHITE
        ).scale(0.9)

        # Agrupar o texto em bloco vertical
        bloco = VGroup(texto_campo_v).arrange(DOWN, aligned_edge=LEFT)

        # Caixa (retângulo) ao redor do texto
        caixa = SurroundingRectangle(
            bloco,
            color=WHITE,
            buff=0.3  # margem interna
        )
        caixa.set_fill(opacity=0.2)

        # Agrupar tudo
        anotacao = VGroup(caixa, bloco)

        self.add_fixed_in_frame_mobjects(anotacao)
        anotacao.to_corner(DOWN)
        self.wait(5)

        self.remove_fixed_in_frame_mobjects(anotacao)
        self.remove(anotacao)

        self.remove(axes, *vetores, reta, plane, dot, v_original, w_original, v(dot), w(dot))
        self.wait()


class TextoCalculoCamposVeWCena01(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Texto na tela, cálculo de campos V e W pi-relacionados com ddx e ddy.
        """

        # Mostrar a expressão da diferencial da projeção pi
        texto_1 = Tex(
            r"A diferencial da função $\pi$ em qualquer ponto $p \in \mathbb{R}^3$, ",
            r"na sua forma matricial, é dada por \\",
            r"$d \pi (p) = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{bmatrix}$",
            color=WHITE
        ).scale(0.7).to_corner(UP)

        texto_2 = Tex(
            r"\begin{flushleft}"
            r"Para um ponto $p$ qualquer, "
            r"\end{flushleft}"
        ).scale(0.7).move_to(UP)

        texto_3 = MathTex(
            r"\begin{aligned}"
            r"V(p)"
            r"&= V_x(p) \frac{\partial}{\partial x} \right| _p +"
            r"V_y(p) \frac{\partial}{\partial y} \right| _p + V_z(p) \frac{\partial}{\partial z} \right| _p\\"
            r"&= V_x(p) \cdot e_1 + V_y(p) \cdot e_2 + V_z(p) \cdot e_3 \\"
            r"&= \begin{bmatrix} V_x(p) \\ V_y(p) \\ V_z(p) \end{bmatrix}"
            r"\end{aligned}"
        ).scale(0.7).move_to(DOWN)

        textos = [texto_1, texto_2, texto_3]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.remove(texto)
            self.wait()
            self.play(Write(texto), run_time=2)
        self.wait(2)

        self.remove(*textos)
        self.wait()

        # Mostrar o cálculo de diferencial \cdot Vp = e_1 e como V = 1 d/dx + 0 d/dy + u d/dz
        texto_1 = Tex(
            r"Nessas condições, podemos escrever"
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"\begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{bmatrix} \cdot"
            r"\begin{bmatrix} V_x(p) \\ V_y(p) \\ V_z(p) \end{bmatrix} = "
            r"\begin{bmatrix} 1 \\ 0 \end{bmatrix}"
        ).scale(0.7).move_to(UP * 2)

        texto_3 = Tex(
            r"\begin{flushleft}"
            r"Deste modo, $V_x \equiv 1$ e $V_y \equiv 0$, enquanto que $V_z$ fica livre para ser uma "
            r"função de classe $C^{\infty}$. Chamemos de $u(x, y, z) = V_z$. Assim, "
            r"\end{flushleft}"
        ).scale(0.7).move_to(UP * 0)

        texto_4 = MathTex(
            r"\begin{aligned}"
            r"V"
            r"&= 1 \cdot \frac{\partial}{\partial x} + 0 \cdot \frac{\partial}{\partial y} +"
            r"u(x,y,z) \frac{\partial}{\partial z} \\"
            r"&= \frac{\partial}{\partial x} + u(x,y,z) \frac{\partial}{\partial z}"
            r"\end{aligned}"
        ).scale(0.7).move_to(DOWN * 2)

        textos = [texto_1, texto_2, texto_3, texto_4]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.remove(texto)
            self.wait()
            self.play(Write(texto), run_time=2)
        self.wait(2)

        self.remove(*textos)
        self.wait()

        texto_1 = Tex(
            r"\begin{flushleft}"
            r"De modo análogo, para encontrarmos um campo $W$ que seja $\pi$-relacionado com "
            r"$\frac{\partial}{\partial y}$ em $\mathbb{R}^2$, devemos achar um campo "
            r"que cumpre\\"
            r"$d \pi (p) \cdot W(p) = \frac{\partial}{\partial y} \right| _p$. "
            r"Semelhante às contas que fizemos para o campo $V$, temos que:"
            r"\end{flushleft}"
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"\begin{aligned}"
            r"W(p)"
            r"&= W_x(p) \frac{\partial}{\partial x} \right| _p +"
            r"W_y(p) \frac{\partial}{\partial y} \right| _p + W_z(p) \frac{\partial}{\partial z} \right| _p\\"
            r"&= W_x(p) \cdot e_1 + W_y(p) \cdot e_2 + W_z(p) \cdot e_3 \\"
            r"&= \begin{bmatrix} W_x(p) \\ W_y(p) \\ W_z(p) \end{bmatrix}"
            r"\end{aligned}"
        ).scale(0.7)

        texto_3 = MathTex(
            r"\begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{bmatrix} \cdot"
            r"\begin{bmatrix} W_x(p) \\ W_y(p) \\ W_z(p) \end{bmatrix} = "
            r"\begin{bmatrix} 0 \\ 1 \end{bmatrix}"
        ).scale(0.7).move_to(DOWN * 3)

        texto_4 = Tex(
            r"\begin{flushleft}"
            r"Deste modo, $W_x \equiv 0$ e $W_y \equiv 1$, enquanto que $W_z$ fica livre para ser uma "
            r"função de classe $C^{\infty}$. Chamemos de $v(x, y, z) = W_z$. Assim, "
            r"\end{flushleft}"
        ).scale(0.7).to_corner(UP)

        texto_5 = MathTex(
            r"\begin{aligned}"
            r"W"
            r"&= 0 \cdot \frac{\partial}{\partial x} + 1 \cdot \frac{\partial}{\partial y} +"
            r"v(x,y,z) \frac{\partial}{\partial z} \\"
            r"&= \frac{\partial}{\partial y} + v(x,y,z) \frac{\partial}{\partial z}"
            r"\end{aligned}"
        ).scale(0.7).move_to(UP)

        textos = [texto_1, texto_2, texto_3, texto_4, texto_5]

        textos_na_tela = []  # Guarda quais textos estão sendo exibidos
        for texto in textos:

            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto, run_time=2))
            self.wait(2)

            textos_na_tela.append(texto)

            if len(textos_na_tela) == 3:  # Máximo de textos por vez
                self.remove(*textos_na_tela)
                textos_na_tela.clear()

        self.remove(*textos_na_tela)
        self.wait(4)


class TextoCalculoCamposVeWCena02(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Texto na tela, cálculo dos campos V e W pi-relacionados com ddx e ddy e subordinados a D.
        """

        # Sabendo as expressões para V e W, fazer continhas para ver que V = -Y e W = X - xY

        texto_1 = Tex(
            r"\begin{flushleft}"
            r"Para que o campo $V = \frac{\partial}{\partial x} + u(x, y, z) \frac{\partial}{\partial z}$ "
            r"seja subordinado à distribuição, devemos ter que $V = fX + gY$ para certas "
            r"$f, g \in C^{\infty}(\mathbb{R}^3)$. Disso, segue que:"
            r"\end{flushleft}"
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"\begin{aligned}"
            r"\frac{\partial}{\partial x} + u(x, y, z) \frac{\partial}{\partial z}"
            r"&= fX + gY\\"
            r"&= f \left(x \frac{\partial}{\partial x} + \frac{\partial}{\partial y} +"
            r"   x(y + 1) \frac{\partial}{\partial z} \right) + "
            r"   g \left(\frac{\partial}{\partial x} + y \frac{\partial}{\partial z} \right)\\"
            r"&= \left( f \cdot x + g \right) \frac{\partial}{\partial x} +"
            r"   \left( f \right) \frac{\partial}{\partial y} +"
            r"   \left( f \cdot x(y + 1) + g \cdot y \right) \frac{\partial}{\partial z}"
            r"\end{aligned}"
        ).scale(0.7)

        texto_3 = Tex(
            r"\begin{center}"
            r"Comparando ambos os lados, segue que \\"
            r"$f \cdot x + g = 1$\\"
            r"$f = 0$\\"
            r"$f \cdot x(y + 1) + g \cdot y = u$\\"
            r"\end{center}"
        ).scale(0.7).move_to(DOWN)

        texto_4 = Tex(
            r"\begin{center}"
            r"Logo, $f \equiv 0$, $g \equiv 1$ e $u(x, y, z) = y$, o que nos permite concluir que\\"
            r"$V = \frac{\partial}{\partial x} + y \frac{\partial}{\partial z} = Y$"
            r"\end{center}"
        ).scale(0.7).move_to(DOWN * 3)

        self.add_fixed_in_frame_mobjects(texto_1)
        self.play(Write(texto_1, run_time=2))
        self.wait(2)
        self.add_fixed_in_frame_mobjects(texto_2)
        self.play(Write(texto_2, run_time=2))
        self.wait(2)

        self.remove(texto_1)
        self.wait()
        self.play(texto_2.animate.to_corner(UP))

        self.add_fixed_in_frame_mobjects(texto_3)
        self.play(Write(texto_3, run_time=2))
        self.wait(2)
        self.add_fixed_in_frame_mobjects(texto_4)
        self.play(Write(texto_4, run_time=2))
        self.wait(8)

        self.remove(texto_1, texto_2, texto_3, texto_4)
        self.wait(2)

        texto_1 = Tex(
            r"\begin{flushleft}"
            r"Analogamente, para que o campo $W = \frac{\partial}{\partial y} +"
            r"v(x, y, z) \frac{\partial}{\partial z}$ "
            r"seja subordinado à distribuição, devemos ter que $W = fX + gY$ para certas "
            r"$f, g \in C^{\infty}(\mathbb{R}^3)$. Disso, segue que:"
            r"\end{flushleft}"
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"\begin{aligned}"
            r"\frac{\partial}{\partial y} + v(x, y, z) \frac{\partial}{\partial z}"
            r"&= fX + gY\\"
            r"&= f \left(x \frac{\partial}{\partial x} + \frac{\partial}{\partial y} +"
            r"   x(y + 1) \frac{\partial}{\partial z} \right) + "
            r"   g \left(\frac{\partial}{\partial x} + y \frac{\partial}{\partial z} \right)\\"
            r"&= \left( f \cdot x + g \right) \frac{\partial}{\partial x} +"
            r"   \left( f \right) \frac{\partial}{\partial y} +"
            r"   \left( f \cdot x(y + 1) + g \cdot y \right) \frac{\partial}{\partial z}"
            r"\end{aligned}"
        ).scale(0.7)

        texto_3 = Tex(
            r"\begin{center}"
            r"Comparando ambos os lados, segue que \\"
            r"$f \cdot x + g = 0$\\"
            r"$f = 1$\\"
            r"$f \cdot x(y + 1) + g \cdot y = v$\\"
            r"\end{center}"
        ).scale(0.7).move_to(DOWN)

        texto_4 = Tex(
            r"\begin{center}"
            r"Logo, $f \equiv 1$, $g = -x$ e $v(x, y, z) = 1 \cdot x(y + 1) - xy = x$, o que nos permite concluir que\\"
            r"$W = \frac{\partial}{\partial y} + x \frac{\partial}{\partial z} = X - xY$"
            r"\end{center}"
        ).scale(0.7).move_to(DOWN * 3)

        self.add_fixed_in_frame_mobjects(texto_1)
        self.play(Write(texto_1, run_time=2))
        self.wait(2)
        self.add_fixed_in_frame_mobjects(texto_2)
        self.play(Write(texto_2, run_time=2))
        self.wait(2)

        self.remove(texto_1)
        self.wait()
        self.play(texto_2.animate.to_corner(UP))

        self.add_fixed_in_frame_mobjects(texto_3)
        self.play(Write(texto_3, run_time=2))
        self.wait(2)
        self.add_fixed_in_frame_mobjects(texto_4)
        self.play(Write(texto_4, run_time=2))
        self.wait(8)

        self.remove(texto_1, texto_2, texto_3, texto_4)
        self.wait(2)

        # Obs: durante a apresentação, falar que, como a projeção restrita à distribuição
        # é isomorfismo, a naturalidade do colchete garante que V e W têm [V, W] = 0

        # Prova real de que [V, W] = 0
        texto_1 = MathTex(
            r"V = \frac{\partial}{\partial x} + y \frac{\partial}{\partial z}"
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"W = \frac{\partial}{\partial y} + x \frac{\partial}{\partial z}"
        ).scale(0.7).move_to(UP * 2)

        texto_3 = MathTex(
            r"\begin{aligned}"
            r"[V, W] "
            r"&= \sum _{j=1} ^3 (V(b_j) - W(a_j)) \frac{\partial}{\partial u_j}\\"
            r"&= (V(0) - W(1)) \frac{\partial}{\partial x} +"
            r"(V(1) - W(0)) \frac{\partial}{\partial y} + "
            r"(V(x) - W(y)) \frac{\partial}{\partial z}\\"
            r"&= (0 - 0) \frac{\partial}{\partial x} +"
            r"(0 - 0) \frac{\partial}{\partial y} + "
            r"(1 - 1) \frac{\partial}{\partial z}\\"
            r"&= 0"
            r"\end{aligned}"
        ).scale(0.7).move_to(DOWN)

        textos = [texto_1, texto_2, texto_3]
        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto, run_time=2))
            self.wait(2)
        self.wait()
        self.remove(*textos)
        self.wait()


class TextoLema02TeoremaFrobenius(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Texto na tela, enunciação do Lema 02 teorema de Frobenius geométrico
        """

        texto_1 = Tex(
            r"\begin{flushleft}"
            r"Lema 2: Sejam $X_1, \dots, X_k$ campos de vetores de classe $C^{\infty}$ "
            r"em um subconjunto aberto $U$ da variedade M. "
            r"Suponha que $\left\{ X_i(q) \right\} _{1 \leq i \leq k}$ seja linearmente independente "
            r"para todo $q \in U$ e que $ \left[ X_i, X_j\right] = 0$ para quaisquer $1 \leq i,j \leq k$.\\"
            r"Então, para todo $p \in U$, existe um sistema de coordenadas "
            r"$\varphi : V \longrightarrow \varphi (V) \subset \mathbb{R}^m$, com $p \in V \subset U$, tal que "
            r"$X_i(q) = \frac{\partial}{\partial u_i} \right| _{q}$ para qualquer $1 \leq i \leq k$."
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).to_corner(UP)

        texto_2 = Tex(
            r"\begin{flushleft}"
            r"O lema 1 nos deu uma forma de obter $k$ campos linearmente independentes que comutam."
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).move_to(DOWN * 0.25)

        texto_3 = Tex(
            r"\begin{flushleft}"
            r"Novamente, a demonstração do lema 2 nos permite construir tal sistema de coordenadas"
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).move_to(DOWN * 1.25)

        textos = [texto_1, texto_2, texto_3]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto, run_time=4))
            self.wait(2)
        self.wait(2)

        for texto in textos:
            self.remove(texto)
        self.wait()


class TextoCalculoFluxosVeW(ThreeDScene):
    def construct(self):

        """
        CENA: Cálculo dos fluxos de V e W
        """

        texto_1 = Tex(
            r"Dado $p = (x_p, y_p, z_p) \in \mathbb{R}^3$, para encontrar o fluxo "
            r"$\phi _p(t) = (x(t), y(t), z(t))$ de $V$ sabemos que, por um lado,"
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"\phi ' _p(t) = (x'(t), y'(t), z'(t))"
        ).scale(0.7).move_to(UP * 2)

        texto_3 = Tex(
            r"Por outro lado, como a imagem do fluxo é uma curva integral, vale que"
        ).scale(0.7).move_to(UP * 1)

        texto_4 = MathTex(
            r"\begin{aligned}"
            r"\phi ' _p(t)"
            r"&= V \left(\phi _p(t) \right)\\"
            r"&= V_x \left(\phi _p(t) \right) \frac{\partial}{\partial x} +"
            r"   V_y \left(\phi _p(t) \right) \frac{\partial}{\partial y} +"
            r"   V_z \left(\phi _p(t) \right) \frac{\partial}{\partial z}\\"
            r"&= \left( V_x (\phi _p(t)), V_y (\phi _p(t)), V_z (\phi _p(t)) \right)\\"
            r"&= \left( V_x (x(t),y(y),z(t)), V_y (x(t),y(y),z(t)), V_z (x(t),y(y),z(t)) \right)"
            r"\end{aligned}"
        ).scale(0.7).move_to(DOWN * 1.5)

        texto_5 = Tex(
            r"Logo, estamos interessados em resolver o seguinte sistema de EDOs"
        ).scale(0.7).to_corner(UP)

        texto_6 = MathTex(
            r"\begin{cases}"
            r"x'(t) = V_x (x(t),y(t),z(t))\\"
            r"y'(t) = V_y (x(t),y(t),z(t))\\"
            r"z'(t) = V_z (x(t),y(t),z(t))\\"
            r"\end{cases}"
        ).scale(0.7).move_to(UP * 2)

        texto_7 = Tex(
            r"Lembrando que $V_x \equiv 1$, $V_y \equiv 0$ e $V_z(x,y,z) = y$, o sistema se torna"
        ).scale(0.7).move_to(UP * 0.5)

        texto_8 = MathTex(
            r"\begin{cases}"
            r"x'(t) = 1\\"
            r"y'(t) = 0\\"
            r"z'(t) = y(t)\\"
            r"\end{cases}"
        ).scale(0.7).move_to(DOWN)

        texto_9 = Tex(
            r"Cuja solução é $x(t) = x_p + t$, $y(t) = y_p$ e $z(t) = z_p + t y_p$. Assim, "
            r"$\phi^V (t, p) = (x + t, y, z + ty)$."
        ).scale(0.7).move_to(DOWN * 3)

        texto_10 = Tex(
            r"Analogamente para o campo $W$, sendo $\phi _p(t) = (x(t), y(t), z(t))$ seu fluxo,"
        ).scale(0.7).to_corner(UP)

        texto_11 = MathTex(
            r"\begin{cases}"
            r"x'(t) = W_x (x(t),y(t),z(t))\\"
            r"y'(t) = W_y (x(t),y(t),z(t))\\"
            r"z'(t) = W_z (x(t),y(t),z(t))\\"
            r"\end{cases}"
        ).scale(0.7).move_to(UP * 2)

        texto_12 = MathTex(
            r"\begin{cases}"
            r"x'(t) = 0\\"
            r"y'(t) = 1\\"
            r"z'(t) = x(t)\\"
            r"\end{cases}"
        ).scale(0.7)

        texto_13 = Tex(
            r"E assim $x(t) = x_p$, $y(t) = y_p + t$ e $z(t) = z_p + t x_p$, de modo que "
            r"$\phi^W (t, p) = (x, y + t, z + tx)$."
        ).scale(0.7).move_to(DOWN * 2)

        textos = [texto_1, texto_2, texto_3, texto_4, texto_5, texto_6,
                  texto_7, texto_8, texto_9, texto_10, texto_11, texto_12, texto_13]

        textos_na_tela = []

        for c in range(4):
            self.add_fixed_in_frame_mobjects(textos[c])
            self.play(Write(textos[c]))
            textos_na_tela.append(textos[c])
            self.wait()
        self.wait(3)
        self.remove(*textos_na_tela)
        textos_na_tela.clear()

        for c in range(4, 9):
            self.add_fixed_in_frame_mobjects(textos[c])
            self.play(Write(textos[c]))
            textos_na_tela.append(textos[c])
            self.wait()
        self.wait(3)
        self.remove(*textos_na_tela)
        textos_na_tela.clear()

        for c in range(9, 13):
            self.add_fixed_in_frame_mobjects(textos[c])
            self.play(Write(textos[c]))
            textos_na_tela.append(textos[c])
            self.wait()
        self.wait(3)
        self.remove(*textos_na_tela)
        textos_na_tela.clear()


class TextoSistemaDeCoordenadas(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        """
        CENA: Ilustração do sistema de coordenadas e das subvariedades integrais
        """

        # Texto para a definição da big_phi (\Phi)
        texto_1 = Tex(
            r"Sabendo que os fluxos de V e W são dados por \\"
            r"$\phi ^V _t (x, y, z) = (x + t, y, z + ty)$ \\"
            r"$\phi ^W _t (x, y, z) = (x, y + t, z + tx)$ \\"
            r"definimos uma função $\Phi : \mathbb{R}^3 \longrightarrow \mathbb{R}^3$ dada por",
            color=WHITE
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"\begin{aligned}"
            r"\Phi(u,v,w)"
            r"&= \phi_u^V \circ \phi_v^W (0,0,w)\\"
            r"&= \phi_u^V (0,v,w)\\"
            r"&= (u,\, v,\, w + uv)"
            r"\end{aligned}",
            color=WHITE
        ).scale(0.7).move_to(UP * 0)

        texto_3 = Tex(
            r"Vejamos um exemplo",
            color=WHITE
        ).scale(0.7).move_to(DOWN * 2)

        textos = [texto_1, texto_2, texto_3]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto), run_time=4)
            self.wait(2)
        self.remove(*textos)
        self.wait(2)


class SistemaDeCoordenadas(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES)

        reta = Line3D(
            start=np.array([0, 0, -10]),
            end=np.array([0, 0, 10]),
            color=WHITE,
            thickness=0.03
        )

        self.add(axes_dom, axes_cod, dot_dom, dot_cod, reta)
        self.add(linha_x_dom, linha_y_dom, linha_z_dom)
        self.add(linha_x_cod, linha_y_cod, linha_z_cod)

        texto_dominio = Tex(
            r"Domínio"
        ).scale(0.8).to_corner(UL)

        texto_contradominio = Tex(
            r"Contradomínio"
        ).scale(0.8).to_corner(UR)

        texto_exemplo = MathTex(
            r"\Phi(1, 2, 3) = \phi_1^V \circ \phi_2^W (0,0,3) = (1,\, 2,\, 3 + 1 \cdot 2) = (1,\, 2,\, 5)",
            color=WHITE
        ).scale(0.8)

        # Agrupar o texto em bloco vertical
        bloco = VGroup(texto_exemplo).arrange(DOWN, aligned_edge=LEFT)

        # Caixa (retângulo) ao redor do texto
        caixa = SurroundingRectangle(
            bloco,
            color=WHITE,
            buff=0.3  # margem interna
        )
        caixa.set_fill(opacity=0.2)

        # Agrupar tudo
        anotacao = VGroup(caixa, bloco)

        self.add_fixed_in_frame_mobjects(texto_dominio)
        self.play(Write(texto_dominio), run_time=2)
        self.wait()
        self.add_fixed_in_frame_mobjects(texto_contradominio)
        self.play(Write(texto_contradominio), run_time=2)
        self.wait(5)

        self.add_fixed_in_frame_mobjects(anotacao)
        anotacao.to_corner(DOWN)
        self.wait(5)

        # Ponto no domínio vagando pelo R^3 (para mostrar as linhas x, y e z no dom e cod)
        self.play(dot_dom.animate.move_to(axes_dom.c2p(0, 0, 3)), run_time=1)
        self.wait(3)
        self.play(dot_dom.animate.move_to(axes_dom.c2p(0, 2, 3)), run_time=1)
        self.wait(3)
        self.play(dot_dom.animate.move_to(axes_dom.c2p(1, 2, 3)), run_time=1)
        self.wait(6)
        self.remove(anotacao)
        self.wait(3)

        self.play(dot_dom.animate.move_to(axes_dom.c2p(2, 2, 0)), run_time=1)
        self.wait()

        # TracedPaths para mostrar o caminho percorrido pelos pontos dot_dom e dot_cod
        traced_dom = TracedPath(
            lambda: dot_dom.get_center(),
            stroke_color=YELLOW,
            stroke_width=4
        )
        traced_cod = TracedPath(
            lambda: dot_cod.get_center(),
            stroke_color=YELLOW,
            stroke_width=4
        )

        # Adiciona os TracedPaths
        self.add(traced_dom)
        axes_cod.add(traced_cod)  # o tracedpath do cod pertence ao sistema de eixos axes_cod para girar com ele

        # Faz o dot_dom se mover em zigue-zague, mantendo sua coordenada z constante
        for (x0, y0, z0) in lista_de_pontos_sistemadecoordenadas:
            self.play(
                dot_dom.animate.move_to(axes_dom.c2p(x0, y0, z0)),
                run_time=0.5,
                rate_func=linear
            )
        self.wait(2)

        self.remove(dot_cod)
        self.wait()

        # Rotaciona o axes_cod para vermos melhor o caminho traçado por dot_cod
        self.play(
            Rotate(
                axes_cod,
                angle=2 * PI,
                axis=OUT,
                run_time=25,
                rate_func=linear
            )
        )

        self.wait(3)

        # Exemplo: $\Phi(1, 2, 3) = \varphi_1^V \circ \varphi_2^W (0, 0, 3)$
        # A partir do ponto (0, 0, 3), andar uma unidade no fluxo de V e duas no fluxo de W


class TextoInverterSistemaDeCoordenadas(ThreeDScene):
    def construct(self):

        # Calcular e inverter $\Phi$, obtendo o sistema de coordenadas. As subvariedades integrais surgem
        # ao igualarmos a terceira componente do sistema de coordenadas a uma constante

        texto_1 = Tex(
            r"Invertendo a função $\Phi$, obtemos"
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"\begin{aligned}"
            r"\Phi (u,v,w)"
            r"&= (u, v, w + uv)\\"
            r"&= (x, y, z)"
            r"\end{aligned}"
        ).scale(0.7).move_to(UP * 2)

        texto_3 = Tex(
            r"Logo, $x = u$, $y = v$ e $z = w + uv$, donde segue que "
        ).scale(0.7).move_to(UP)

        texto_4 = MathTex(
            r"\Phi ^{-1} (x,y,z) = (u,v,w) = (x,y,z - xy)"
        ).scale(0.7)

        texto_5 = Tex(
            r"Assim, as subvariedades integrais da distribuição são curvas de nível do tipo $z - xy = C$."
        ).scale(0.7).move_to(DOWN)

        textos = [texto_1, texto_2, texto_3, texto_4, texto_5]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto))
            self.wait(2)
        self.wait(2)
        self.remove(*textos)
        self.wait()


class GraficoComCorte(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False
        self.set_camera_orientation(phi=60*DEGREES, theta=45*DEGREES)

        self.add(axes)
        self.wait()

        def superf(u, v):
            x = u
            y = v
            z = x**2/2 - y**2/2
            return np.array([x, y, z])

        superficie = Surface(
            superf,
            u_range=[-3, 3],
            v_range=[-3, 3],
            resolution=(20, 20),
            fill_opacity=0.7,
            fill_color=BLUE,
            stroke_color=WHITE
        )

        self.add(superficie)
        self.wait(1)

        def curva(t):
            x = t
            y = 1
            z = x ** 2 / 2 - 1 ** 2 / 2
            return np.array([x, y, z])

        interseccao = ParametricFunction(
            curva,
            t_range=(-3, 3),
            color=YELLOW,
            stroke_width=6
        )

        self.play(Create(interseccao))
        self.wait(2)

        # Superfície cortada (apenas y >= 1)
        superficie_cortada = Surface(
            superf,
            u_range=[-3, 3],
            v_range=[1, 3],  # << CORTE AQUI
            resolution=(10, 10),
            fill_opacity=0.7,
            fill_color=BLUE,
            stroke_color=WHITE
        )

        # Remove a superfície inteira e coloca apenas a parte "para frente" (y>=1)
        self.remove(superficie)
        self.add(superficie_cortada)
        self.wait(2)

        self.move_camera(
            phi=90 * DEGREES,
            theta=-90 * DEGREES,
            run_time=1,
            rate_func=linear
        )
        self.wait()

        p = np.array([1, 1, 0])
        dot_p = Dot3D(p, color=WHITE)
        self.add(dot_p)

        vetor_dx = Arrow3D(
            start=p,
            end=p + np.array([1, 0, 0]),
            color=RED,
            thickness=0.04
        )
        self.play(Create(vetor_dx))
        self.wait(1)

        texto_p = MathTex(r"p", color=WHITE).scale(0.7)
        texto_dx = MathTex(r"\frac{\partial}{\partial x}\Big|_p", color=WHITE).scale(0.6)

        # Posicionamento
        texto_p.move_to(RIGHT + UP * 0.4)
        texto_dx.move_to(RIGHT * 1.5 + DOWN * 0.5)

        self.add_fixed_in_frame_mobjects(texto_p, texto_dx)
        self.play(Write(texto_p), Write(texto_dx))
        self.wait(1)

        gamma_prime = np.array([1, 0, 1])  # derivada da curva em x=1
        reta_tg = Line3D(
            start=p - 2 * gamma_prime,
            end=p + 2 * gamma_prime,
            color=ORANGE,
            thickness=0.02
        )

        self.play(Create(reta_tg))
        self.wait(1)

        origem2 = np.array([2, 1, 0])
        ponto_tg = p + gamma_prime  # (2,1,1)
        vetor_vertical = Arrow3D(
            start=origem2,
            end=ponto_tg,
            color=BLUE,
            thickness=0.04
        )

        self.play(Create(vetor_vertical))
        self.wait(1)

        texto_du = MathTex(
            r"\frac{\partial u}{\partial x}(p)\cdot\frac{\partial}{\partial z}",
            color=WHITE
        ).scale(0.6)

        texto_du.move_to(RIGHT * 2.8 + UP * 0.5)

        self.add_fixed_in_frame_mobjects(texto_du)
        self.play(Write(texto_du))

        self.wait(3)

        self.remove(axes, interseccao, superficie_cortada, dot_p, reta_tg,
                    vetor_dx, vetor_vertical, texto_p, texto_dx, texto_du)

        self.wait(2)

        texto_1 = Tex(
            r"\begin{flushleft}"
            r"Sendo $u : U \subset \mathbb{R}^2 \longrightarrow \mathbb{R}$ de classe $C^{\infty}$ e definindo "
            r"$V = \frac{\partial}{\partial x} + \frac{\partial u}{\partial x} \frac{\partial}{\partial z}$, "
            r"sabemos que $V(p) \in T_p Gr(u)$ para todo $p \in Gr(u)$."
            r"\end{flushleft}"
        ).scale(0.7).to_corner(UP)

        texto_2 = Tex(
            r"\begin{flushleft}"
            r"De modo análogo, definindo "
            r"$W = \frac{\partial}{\partial y} + \frac{\partial u}{\partial y} \frac{\partial}{\partial z}$, "
            r"temos que $W(p) \in T_p Gr(u)$ para todo $p \in Gr(u)$."
            r"\end{flushleft}"
        ).scale(0.7).move_to(UP)

        texto_3 = Tex(
            r"\begin{flushleft}"
            r"Assim, a distribuição gerada por $V$ e $W$ é sempre tangente ao gráfico de $u$."
            r"\end{flushleft}"
        ).scale(0.7).move_to(DOWN)

        textos = [texto_1, texto_2, texto_3]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto, run_time=2))
            self.wait(2)
        self.wait()
        self.remove(*textos)


class TextoSistemaDeEDPs(ThreeDScene):
    def construct(self):

        texto_1 = Tex(
            r"Considere o seguinte sistema de EDPs"
        ).scale(0.7).to_corner(UP)

        texto_2 = MathTex(
            r"\begin{cases}"
            r"\frac{\partial u}{\partial x}(x,y) = \alpha (x, y, u(x,y))\\"
            r"\frac{\partial u}{\partial y}(x,y) = \beta (x, y, u(x,y))\\"
            r"u(x_0, y_0) = z_0"
            r"\end{cases}"
        ).scale(0.7).move_to(UP * 2)

        texto_3 = Tex(
            r"\begin{center}"
            r"Uma condição necessária para que esse sistema admita solução é: \\"
            r"$\frac{\partial ^2 u}{\partial x \partial y} = \frac{\partial ^2 u}{\partial y \partial x}$\\"
            r"Isto é, que \\"
            r"$\frac{\partial \alpha}{\partial y} + \beta \frac{\partial \alpha}{\partial z} = "
            r"\frac{\partial \beta}{\partial x} + \alpha \frac{\partial \beta}{\partial z}$"
            r"\end{center}"
        ).scale(0.7).move_to(DOWN)

        texto_4 = Tex(
            r"A conexão entre o teorema de Frobenius geométrico e a resolução desse sistema é que, caso "
            r"$\alpha$ e $\beta$ satisfaçam a condição acima, os campos \\"
            r"$V = \frac{\partial}{\partial x} + \alpha \frac{\partial}{\partial z} \approx "
            r"\frac{\partial}{\partial x} + \frac{\partial u}{\partial x} \frac{\partial}{\partial z}$\\"
            r"$W = \frac{\partial}{\partial y} + \beta \frac{\partial}{\partial z} \approx "
            r"\frac{\partial}{\partial y} + \frac{\partial u}{\partial y} \frac{\partial}{\partial z}$\\"
            r"irão sempre gerar uma distribuição involutiva (que será sempre tangente ao gráfico de $u(x, y, z)$)."
        ).scale(0.7)

        texto_5 = Tex(
            r"\begin{flushleft}"
            r"Como toda distribuição involutiva é integrável, pelo ponto $(x_0, y_0, z_0)$ passará uma "
            r"subvariedade integral que irá coincidir com o gráfico de $u(x, y, z)$. Logo, o sistema "
            r"admitirá uma (única) solução ($u$)."
            r"\end{flushleft}"
        ).scale(0.7)

        texto_6 = Tex(
            r"No exemplo trabalhado,"
        ).scale(0.7).to_corner(UP)

        texto_7 = MathTex(
            r"V = \frac{\partial}{\partial x} + y \frac{\partial}{\partial z}"
        ).scale(0.7).move_to(UP * 2)

        texto_8 = MathTex(
            r"W = \frac{\partial}{\partial y} + x \frac{\partial}{\partial z}"
        ).scale(0.7).move_to(UP)

        texto_9 = MathTex(
            r"\begin{cases}"
            r"\frac{\partial u}{\partial x}(x,y) = y\\"
            r"\frac{\partial u}{\partial y}(x,y) = x\\"
            r"u(x_0, y_0) = z_0"
            r"\end{cases}"
        ).scale(0.7).move_to(DOWN)

        texto_10 = Tex(
            r"Única solução: $u(x,y) = xy + (z_0 - x_0 y_0)$"
        ).scale(0.7).move_to(DOWN * 3)

        self.add_fixed_in_frame_mobjects(texto_1)
        self.play(Write(texto_1), run_time=1)
        self.wait(2)
        self.add_fixed_in_frame_mobjects(texto_2)
        self.play(Write(texto_2), run_time=2)
        self.wait(2)
        self.add_fixed_in_frame_mobjects(texto_3)
        self.play(Write(texto_3), run_time=2)
        self.wait(4)

        self.remove(texto_1, texto_2)
        self.wait()
        self.play(texto_3.animate.to_corner(UP))
        self.wait()

        self.add_fixed_in_frame_mobjects(texto_4)
        self.play(Write(texto_4), run_time=2)
        self.wait(4)

        self.remove(texto_3)
        self.wait()
        self.play(texto_4.animate.to_corner(UP))
        self.wait()

        self.add_fixed_in_frame_mobjects(texto_5)
        self.play(Write(texto_5), run_time=2)
        self.wait(4)
        self.remove(texto_4, texto_5)
        self.wait()

        textos = [texto_6, texto_7, texto_8, texto_9, texto_10]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.play(Write(texto), run_time=2)
            self.wait(2)
        self.wait(3)
        self.remove(*textos)
