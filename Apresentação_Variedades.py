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


def subvariedade_integral(point, x_range=(-6, 6), y_range=(-6, 6)):
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
        resolution=(30, 30),
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


# ----------------------------------------------
# Criação dos elementos utilizados nas animações
# ----------------------------------------------


# Axes (eixos) principais a serem utilizados em todas as animações
axes = ThreeDAxes(
    x_range=[-6, 6, 1],
    y_range=[-6, 6, 1],
    z_range=[-6, 6, 1],
    x_length=12,
    y_length=12,
    z_length=12,
)

# Ponto que será constantemente utilizado
dot = Dot3D(np.array([-1, -1, -1]), radius=0.08, color=YELLOW)

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
                resolution=2,
                stroke_width=1,
            )
            arrowY = Arrow3D(
                start=p,
                end=p + 0.5 * Yv / norm_Yv,
                color=GREEN,
                resolution=2,
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


class Animation(ThreeDScene):
    def construct(self):

        self.renderer.camera.shading = False

        # Texto do campo X
        texto_campo_X = MathTex(
            r"X = x\,\frac{\partial}{\partial x} + \frac{\partial}{\partial y} + x(y+1)\,\frac{\partial}{\partial z}",
            color=RED
        ).scale(0.7).to_corner(UL)

        # Texto do campo Y
        texto_campo_Y = MathTex(
            r"Y = \frac{\partial}{\partial x} + y\,\frac{\partial}{\partial z}",
            color=GREEN
        ).scale(0.7).to_corner(UL)

        """
        Criação dos elementos:
        """



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

        self.remove(dot, plane)
        self.wait(3)

        """
        CENA: Ilustração dos dois campos de vetores X e Y que geram a distribuição
        """

        # Mostra o campo de vetores X
        self.add(*vector_field_X)
        self.wait(2)
        self.add_fixed_in_frame_mobjects(texto_campo_X)
        self.add(texto_campo_X)
        self.wait(2)
        self.move_camera(
            phi=60 * DEGREES,
            theta=405 * DEGREES,
            run_time=25,
            rate_func=linear
        )
        self.wait()
        self.remove(*vector_field_X, texto_campo_X)
        self.wait()

        # Mostra o campo de vetores Y
        self.add(*vector_field_Y)
        self.wait(2)
        self.add_fixed_in_frame_mobjects(texto_campo_Y)
        self.add(texto_campo_Y)
        self.wait(2)
        self.move_camera(
            phi=60 * DEGREES,
            theta=405 * DEGREES,
            run_time=25,
            rate_func=linear
        )
        self.wait()
        self.remove(*vector_field_Y, texto_campo_Y)
        self.wait()

        """
        CENA: Mostra como os campos geram a distribuição
        """

        # Mostra os campos onde o ponto passa
        self.add(*mini_vector_field_X, *mini_vector_field_Y)

        # --- Movement in straight segments ---
        self.add(dot, plane)
        self.wait()
        for ponto in lista_de_pontos:
            self.play(dot.animate.move_to(ponto), run_time=2, rate_func=smooth)
            self.wait()
        self.wait(5)
        self.remove(axes, dot, plane, *mini_vector_field_Y, *mini_vector_field_Y)
        self.wait()

        """
        CENA: Texto na tela, demonstração do lema para distribuição involutiva
        """

        texto_lema1_1 = Tex(
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
            r"$V(p), W(p) \in D(p)$ para todo $p \in M$. Mostremos que $\left[ X, Y \right](p) \in D(p)$, "
            r"para todo ponto $p \in M$."
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).move_to(UP)

        dem2 = Tex(
            r"\begin{flushleft}"
            r"Sabemos que $V = f_1 \cdot X + g_1 \cdot Y$ e $W = f_2 \cdot X + g_2 \cdot Y$, "
            r"para certas $f_1, f_2, g_1, g_2 \in C^{\infty} (M)$. Desta forma, sendo $p \in M$ qualquer, temos que "
            r"\end{flushleft}"
        ).scale(0.7).move_to(UP)

        dem3 = Tex(
            r"\text{Usando a bilinearidade do colchete de Lie, temos}",
            color=WHITE
        ).scale(0.7).move_to(UP)

        dem4 = MathTex(
            r"\begin{aligned}"
            r"[V, W]"
            r"&= [f_1 X + g_1 Y, f_2 X + g_2 Y] \\"
            r"&= [f_1 X, f_2 X] + [f_1 X, g_2 Y] \\"
            r"& + [g_1 Y, f_2 X] + [g_1 Y, g_2 Y]"
            r"\end{aligned}",
            color=WHITE
        ).scale(0.7).move_to(UP)

        dem5 = Tex(
            r"Pela identidade $[fX, gY] = fg[X,Y] + f(Xg)Y - g(Yf)X$, obtemos:",
            color=WHITE
        ).scale(0.7).move_to(UP)

        dem6 = MathTex(
            r"[f_1 X, f_2 X] = f_1 f_2 [X,X] + f_1(X f_2)X - f_2(X f_1)X, \\[6pt]"
            r"[f_1 X, g_2 Y] = f_1 g_2 [X,Y] + f_1(X g_2)Y - g_2(Y f_1)X, \\[6pt]"
            r"[g_1 Y, f_2 X] = g_1 f_2 [Y,X] + g_1(Y f_2)X - f_2(X g_1)Y, \\[6pt]"
            r"[g_1 Y, g_2 Y] = g_1 g_2 [Y,Y] + g_1(Y g_2)Y - g_2(Y g_1)Y.",
            color=WHITE
        ).scale(0.7).move_to(UP)

        dem7 = Tex(
            r"Assim, temos $[V, W] = \alpha \cdot X + \beta \cdot Y + \theta \cdot [X, Y]$, com ",
            color=WHITE
        ).scale(0.7).move_to(UP)

        dem8 = Tex(
            r"$\alpha = f_1 X(f_2) - f_2 X(f_1) - g_2 Y(f_1) + g_1 Y(f_2)$;\\[6pt]"
            r"$\beta  = f_1 X(g_2) - f_2 X(g_1) + g_1 Y(f_2) - g_2 Y(g_1)$;\\[6pt]"
            r"$\theta = f_1 g_2 - g_1 f_2$.",
            color=WHITE
        ).scale(0.7).move_to(UP)

        dem9 = Tex(
            r"Logo, para todo $p \in M$, "
            r"$[V, W](p) = \alpha(p)X(p) + \beta(p)Y(p) + \theta(p)[X, Y](p) \in D(p)$. \\"
            r"Portanto, D é involutiva.",
            color=WHITE
        ).scale(0.7).move_to(UP)

        textos = [texto_lema1_1, dem1, dem2, dem3, dem4, dem5, dem6, dem7, dem8, dem9]

        for texto in textos:
            self.add_fixed_in_frame_mobjects(texto)
            self.remove(texto)

        for texto in textos:
            self.play(Write(texto, run_time=2))
            self.wait(2)
            self.remove(texto)
        self.wait()

        """
        CENA: Texto na tela, cálculo do colchete de X e Y
        """

        # OBSERVAÇÃO: FAZER O CÁLCULO

        """
        CENA: Texto na tela, enunciação do teorema geométrico de Frobenius
        """

        # OBSERVAÇÃO: ENUNCIAR

        """
        CENA: Mostrar subvariedades integráveis como curvas de nível
        """

        self.wait()
        self.add(dot, plane)
        self.wait()
        self.move_camera(
            phi=60 * DEGREES,
            theta=120 * DEGREES,
            run_time=1,
            rate_func=smooth
        )
        self.wait()



        self.add(curvas_nivel_subv_integraveis)


        for curve in [level_path_01, level_path_02, level_path_03, path_04]:
            self.play(MoveAlongPath(dot, curve), run_time=4, rate_func=smooth)
            self.wait(2)
        self.wait(4)

        foliacoes = []

        for z in [-6, 6, 12]:
            point = np.array([1, 2, z])
            foliacoes.append(
                subvariedade_integral(point)
            )

        self.remove(plane)

        for folha in foliacoes:
            self.add(folha)
            self.wait(1)
        self.wait(3)

        self.remove(dot, curvas_nivel_subv_integraveis)
        for folha in foliacoes:
            self.remove(folha)
        self.wait(5)

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

        for texto in cena_lema_frobenius_1:
            self.play(Write(texto, run_time=2))
            self.wait(2)
        self.wait(2)

        for texto in cena_lema_frobenius_1:
            self.remove(texto)
        self.wait()

        """
        CENA: Campos V e W pi-relacionados com d/dx e d/dy
        """

        self.wait()

        self.add(dot, plane, axes)
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
            run_time=2,
            rate_func=linear
        )
        self.wait()

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
                resolution=2,
                stroke_width=1,
            )

            vetores.append(vetor)

        for vetor in vetores:
            self.add(vetor)
            self.wait()
        self.wait(2)

        reta = Line3D(
            start=np.array([2, 1, -6]),
            end=np.array([2, 1, 6]),
            color=YELLOW,
            thickness=0.01
        )
        self.add(reta)
        self.wait(6)

        self.remove(axes, *vetores, reta, plane, dot, v_original, w_original, v(dot), w(dot))

        # $d \pi (p) \cdot V_p = \frac{\partial}{\partial x}|_p = e_1$, for all $p \in M$

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
        self.wait(2)

        self.remove_fixed_in_frame_mobjects(anotacao)
        self.remove(anotacao)

        # Mostrar a expressão da diferencial da projeção pi

        texto_1 = Tex(
            r"A diferencial da função $\pi$ em qualquer ponto $p$, ",
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
            r"V_y(p) \frac{\partial}{\partial y} \right| _p + V_z(p) \frac{\partial}{\partial z} \\"
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

        self.remove(texto_1, texto_2, texto_3, texto_4)
        self.wait()

        # Analogamente, mostrar essas animações para w_p

        # Sabendo as expressões para V e W, fazer continhas para ver que V = -Y e W = X - xY

        # Falar que, como a projeção restrita à distribuição é iso, a naturalidade do colchete
        # garante que V e W tem [V, W] = 0

        """
        CENA: Texto na tela, enunciação do Lema 02 teorema de Frobenius geométrico
        """

        frobenius_lema2_1 = Tex(
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

        frobenius_lema2_2 = Tex(
            r"\begin{flushleft}"
            r"O lema 1 nos deu uma forma de obter $k$ campos linearmente independentes que comutam."
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).move_to(DOWN * 0.25)

        frobenius_lema2_3 = Tex(
            r"\begin{flushleft}"
            r"Novamente, a demonstração do lema 2 nos permite construir tal sistema de coordenadas"
            r"\end{flushleft}",
            color=WHITE
        ).scale(0.7).move_to(DOWN * 1.25)

        cena_lema_frobenius_2 = [frobenius_lema2_1, frobenius_lema2_2, frobenius_lema2_3]

        for texto in cena_lema_frobenius_2:
            self.add_fixed_in_frame_mobjects(texto)

        for texto in cena_lema_frobenius_2:
            self.play(Write(texto))
            self.wait(2)
        self.wait(2)

        for texto in cena_lema_frobenius_2:
            self.remove(texto)
        self.wait()

        """
        CENA: Cálculo dos fluxos de V e W
        """

        # Fazer uma animação do fluxo de V que sai do ponto (1, 1, 1) e vai até (3, 1, 3)

        # No caminho, mostrar o campo V nos pontos (1, 1, 1), (2, 1, 2) e (3, 1, 3) (que serão tangentes ao fluxo)

        # $\phi_p (t) = (x(t), y(t), z(t))$
        # Por um lado, $\phi ' _p (t) = (x'(t), y'(t), z'(t))$
        # Por outro, $\phi ' _p (t) = V(\phi_p(t)) = a_1(\phi_p(t)) \frac{\partial}{\partial x} +
        # a_2(\phi_p(t)) \frac{partial}{\partial y} + a_3(\phi_p(t)) \frac{\partial}{\partial z}$

        # Daí, surge o sistema de equações

        # Mostrar resoluções dos sistemas, mas rapidinho

        # Analogamente para o fluxo de W

        """
        CENA: Ilustração do sistema de coordenadas e das subvariedades integrais
        """

        # Criação dos elementos da cena

        # Texto para a definição da big_phi
        texto_1 = Tex(
            r"Sabendo que os fluxos de V e W são dados por"
            r"$$$$",
            color=WHITE
        ).scale(0.7).to_corner(UP)

        # Cria dois eixos coordenados
        axes_dom = ThreeDAxes(  # axes_dom é o domínio da função big_phi
            x_range=[-5, 5, 1],
            y_range=[-5, 5, 1],
            z_range=[-5, 5, 1],
            x_length=9,
            y_length=9,
            z_length=9,
        ).shift(RIGHT * 3 + DOWN * 3)

        axes_cod = ThreeDAxes(  # axes_cod é o contradomínio da função big_phi
            x_range=[-5, 5, 1],
            y_range=[-5, 5, 1],
            z_range=[-5, 5, 1],
            x_length=9,
            y_length=9,
            z_length=9,
        ).shift(LEFT * 3 + UP * 3)

        # Cria ponto dot_dom (no domínio da big_phi) e sua imagem, dot_cod
        dot_dom = Dot3D(axes_dom.c2p(0, 0, 0), color=YELLOW)  # dom = domínio
        dot_cod = Dot3D(axes_cod.c2p(0, 0, 0), color=YELLOW)  # cod = contradomínio

        # updater do ponto da imagem (dot_cod)
        def update_cod(ponto):
            # pegar coordenadas em eixos do domínio
            u, v, w = axes_dom.p2c(dot_dom.get_center())
            x, y, z = big_phi(u, v, w)
            ponto.move_to(axes_cod.c2p(x, y, z))

        dot_cod.add_updater(update_cod)

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

        # Início da cena

        # Definir $\Phi(u, v, w) = \varphi_u^V \circ \varphi_v^W (0, 0, w)$

        self.add(axes_dom, axes_cod, dot_dom, dot_cod)
        self.add(linha_x_dom, linha_y_dom, linha_z_dom)
        self.add(linha_x_cod, linha_y_cod, linha_z_cod)

        # Ponto no domínio vagando pelo R^3 (para mostrar as linhas x, y e z no dom e cod)
        self.play(dot_dom.animate.move_to(axes_dom.c2p(0, 0, 3)), run_time=2)
        self.play(dot_dom.animate.move_to(axes_dom.c2p(0, 2, 3)), run_time=2)
        self.play(dot_dom.animate.move_to(axes_dom.c2p(1, 2, 3)), run_time=2)
        self.play(dot_dom.animate.move_to(axes_dom.c2p(2, 2, 0)), run_time=1)

        self.wait(5)

        pontos = [
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

        self.add(traced_dom)

        axes_cod.add(traced_cod)  # o tracedpath do cod pertence ao sistema de eixos axes_cod para girar com ele

        for (x, y, z) in pontos:
            self.play(
                dot_dom.animate.move_to(axes_dom.c2p(x, y, z)),
                run_time=0.1,
                rate_func=linear
            )

        self.wait(2)

        self.remove(dot_cod)
        self.play(
            Rotate(
                axes_cod,
                angle=2 * PI,
                axis=OUT,
                run_time=4,
                rate_func=linear
            )
        )

        # Definir $\Phi(u, v, w) = \varphi_u^V \circ \varphi_v^W (0, 0, w)$

        # Exemplo: $\Phi(1, 2, 3) = \varphi_1^V \circ \varphi_2^W (0, 0, 3)$
        # A partir do ponto (0, 0, 3), andar uma unidade no fluxo de V e duas no fluxo de W

        # Mostrar animação com dois R^3's e a aplicação $\Phi$ agindo entre eles (com as cores)

        # Animação fixando w no domínio e variando u e v, obtendo assim as subvariedades integrais

        # Calcular e inverter $\Phi$, obtendo o sistema de coordenadas. As curvas integrais surgem
        # ao igualarmos a terceira componente do sistema de coordenadas a uma constante
