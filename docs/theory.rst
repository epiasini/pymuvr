Mathematical methods
====================

We define a compact formalism for multiunit spike-train metrics by
opportunely characterising the space of multiunit feature vectors as a
tensor product. Previous results from Houghton and Kreuz (2012, *On the
efficient calculation of van Rossum distances*. Network: Computation in
Neural Systems, 2012, 23, 48-58) on a clever formula for multiunit Van
Rossum metrics are then re-derived within this framework, also fixing
some errors in the original calculations.

A compact formalism for kernel-based multiunit spike train metrics
------------------------------------------------------------------

Consider a network with :math:`C` cells. Let

.. math:: \mathcal{U}= \left\{ \boldsymbol{u}^1, \boldsymbol{u}^2, \ldots, \boldsymbol{u}^C \right \}

be an *observation of network activity*, where

.. math:: \boldsymbol{u}^i = \left\{ u_1^i, u_2^i, \ldots, u_{N_{\boldsymbol{u}^i}}^i \right \}

is the (ordered) set of times of the spikes emitted by cell :math:`i`.
Let
:math:`\mathcal{V}= \left\{\boldsymbol{v}^1, \boldsymbol{v}^2, \ldots, \boldsymbol{v}^C\right\}`
be another observation, different in general from :math:`\mathcal{U}`.

To compute a kernel based multiunit distance between :math:`\mathcal{U}`
and :math:`\mathcal{V}`, we map them to the tensor product space
:math:`\mathcal{S} \doteq
\mathbb{R}^C\bigotimes L_2(\mathbb{R}\rightarrow\mathbb{R})` by defining

.. math:: |\mathcal{U}\rangle = \sum_{i=1}^C |i\rangle \otimes |\boldsymbol{u}^i\rangle

where we consider :math:`\mathbb{R}^C` and
:math:`L_2(\mathbb{R}\rightarrow\mathbb{R})` to be equipped with the
usual euclidean distances, consequently inducing an euclidean metric
structure on :math:`\mathcal{S}` too.

Conceptually, the set of vectors
:math:`\left\{ |i\rangle \right\}_{i=1}^C
\subset \mathbb{R}^C` represents the different cells, while each
:math:`|\boldsymbol{u}^i\rangle \in L_2(\mathbb{R}\rightarrow\mathbb{R})`
represents the convolution of a spike train of cell :math:`i` with a
real-valued feature function
:math:`\phi: \mathbb{R}\rightarrow\mathbb{R}`,

.. math:: \langle t|\boldsymbol{u}\rangle = \sum_{n=1}^{N_{\boldsymbol{u}}}\phi(t-u_n)

In practice, we will never use the feature functions directly, but we
will be only interested in the inner products of the :math:`|i\rangle`
and :math:`|\boldsymbol{u}\rangle` vectors. We call
:math:`c_{ij}\doteq\langle i|j \rangle=\langle i|j \rangle_{\mathbb{R}^C}=c_{ji}`
the *multiunit mixing coefficient* for cells :math:`i` and :math:`j`,
and
:math:`\langle \boldsymbol{u}|\boldsymbol{v}\rangle=\langle \boldsymbol{u}|\boldsymbol{v}\rangle_{L_2}`
the *single-unit inner product*,

.. math::

   \label{eq:singleunit_intprod}
     \begin{split}
       \langle \boldsymbol{u}|\boldsymbol{v}\rangle & = \langle \left\{ u_1, u_2, \ldots,
           u_{N}\right\}|\left\{ v_1, v_2, \ldots, v_{M}\right\} \rangle = \\
       &= \int\textrm{d }\!t \langle \boldsymbol{u}|t \rangle\langle t|\boldsymbol{v}\rangle = \int\textrm{d }\!t
       \sum_{n=1}^N\sum_{m=1}^M\phi\left(t-u_n\right)\phi\left(t-v_m\right)\\
       &\doteq \sum_{n=1}^N\sum_{m=1}^M \mathcal{K}(u_n,v_m)
     \end{split}

where :math:`\mathcal{K}(t_1,t_2)\doteq\int\textrm{d }\!t
\left[\phi\left(t-t_1\right)\phi\left(t-t_2\right)\right]` is the
*single-unit metric kernel*, and where we have used the fact that the
feature function :math:`\phi` is real-valued. It follows immediately
from the definition above that
:math:`\langle \boldsymbol{u}|\boldsymbol{v}\rangle=\langle \boldsymbol{v}|\boldsymbol{u}\rangle`.

Note that, given a cell pair :math:`(i,j)` or a spike train pair
:math:`(\boldsymbol{u},\boldsymbol{v})`, :math:`c_{ij}` does not depend
on spike times and :math:`\langle \boldsymbol{u}|\boldsymbol{v}\rangle`
does not depend on cell labeling.

With this notation, we can define the *multi-unit spike train distance*
as

.. math::

   \label{eq:distance_as_intprod}
     \left\Vert |\mathcal{U}\rangle - |\mathcal{V}\rangle \right\Vert^2 = \langle \mathcal{U}|\mathcal{U}\rangle + \langle \mathcal{V}|\mathcal{V}\rangle
   - 2 \langle \mathcal{U}|\mathcal{V}\rangle

where the *multi-unit spike train inner product*
:math:`\langle \mathcal{V}|\mathcal{U}\rangle` between
:math:`\mathcal{U}` and :math:`\mathcal{V}` is just the natural bilinear
operation induced on :math:`\mathcal{S}` by the tensor product
structure:

.. math::

   \label{eq:multiunit_intprod}
     \begin{split}
       \langle \mathcal{V}|\mathcal{U}\rangle &= \sum_{i,j=1}^C
       \langle i|j \rangle\langle \boldsymbol{v}^i|\boldsymbol{u}^j \rangle
       = \sum_{i,j=1}^C c_{ij}\langle \boldsymbol{v}^i|\boldsymbol{u}^i \rangle\\
       &= \sum_{i=1}^C\left[ c_{ii}
         \langle \boldsymbol{v}^i|\boldsymbol{u}^i \rangle + c_{ij}
         \left(\sum_{j<i}\langle \boldsymbol{v}^i|\boldsymbol{u}^j \rangle +
           \sum_{j>i}\langle \boldsymbol{v}^i|\boldsymbol{u}^j \rangle \right) \right]
     \end{split}

But :math:`c_{ij}=c_{ji}` and
:math:`\langle \boldsymbol{v}|\boldsymbol{u}\rangle=\langle \boldsymbol{u}|\boldsymbol{v}\rangle`,
so

.. math::

   \begin{split}
       \sum_{i=1}^C\sum_{j<i}c_{ij}\langle \boldsymbol{v}^i|\boldsymbol{u}^j \rangle &=
       \sum_{j=1}^C\sum_{i<j}c_{ji}\langle \boldsymbol{v}^j|\boldsymbol{u}^i \rangle = 
       \sum_{i=1}^C\sum_{j>i}c_{ji}\langle \boldsymbol{v}^j|\boldsymbol{u}^i \rangle =
       \sum_{i=1}^C\sum_{j>i}c_{ij}\langle \boldsymbol{v}^j|\boldsymbol{u}^i \rangle\\
       &=\sum_{i=1}^C\sum_{j>i}c_{ij}\langle \boldsymbol{u}^i|\boldsymbol{v}^j \rangle
     \end{split}

and

.. math::

   \label{eq:multiprod_j_ge_i}
       \langle \mathcal{V}|\mathcal{U}\rangle = \sum_{i=1}^C\left[ c_{ii}
         \langle \boldsymbol{v}^i|\boldsymbol{u}^i \rangle + c_{ij}
         \sum_{j>i}\left(\langle \boldsymbol{v}^i|\boldsymbol{u}^j \rangle +
           \langle \boldsymbol{u}^i|\boldsymbol{v}^j \rangle \right) \right]

Now, normally we are interested in the particular case where
:math:`c_{ij}` is the same for all pair of distinct cells:

.. math::

   c_{ij} = \begin{cases}
     1 & \textrm{if } i=j\\
     c & \textrm{if } i\neq j
   \end{cases}

and under this assumption we can write

.. math::

   \label{eq:multiprod_constant_c}
     \langle \mathcal{V}|\mathcal{U}\rangle = \sum_{i=1}^C\left[\langle \boldsymbol{v}^i|\boldsymbol{u}^i \rangle + c\sum_{j>i}\left(\langle \boldsymbol{v}^i|\boldsymbol{u}^j \rangle + \langle \boldsymbol{u}^i|\boldsymbol{v}^j \rangle \right) \right]

and

.. math::

   \begin{split}
       \left\Vert |\mathcal{U}\rangle - |\mathcal{V}\rangle \right\Vert^2 =
       \sum_{i=1}^C\Bigg\{&\langle \boldsymbol{u}^i|\boldsymbol{u}^i \rangle +
       c\sum_{j>i}\left(\langle \boldsymbol{u}^i|\boldsymbol{u}^j \rangle + \langle \boldsymbol{u}^i|\boldsymbol{u}^j \rangle
       \right) +\\
       +&\langle \boldsymbol{v}^i|\boldsymbol{v}^i \rangle + c\sum_{j>i}\left(\langle \boldsymbol{v}^i|\boldsymbol{v}^j \rangle +
         \langle \boldsymbol{v}^i|\boldsymbol{v}^j \rangle \right) + \\
       -2&\left[\langle \boldsymbol{v}^i|\boldsymbol{u}^i \rangle +
         c\sum_{j>i}\left(\langle \boldsymbol{v}^i|\boldsymbol{u}^j \rangle + \langle \boldsymbol{u}^i|\boldsymbol{v}^j \rangle
         \right)\right] \Bigg\}
     \end{split}

Rearranging the terms

.. math::

   \begin{gathered}
     \label{eq:multidist_constant_c}
     \left\Vert |\mathcal{U}\rangle - |\mathcal{V}\rangle \right\Vert^2 =
     \sum_{i=1}^C\Bigg[\langle \boldsymbol{u}^i|\boldsymbol{u}^i \rangle + \langle \boldsymbol{v}^i|\boldsymbol{v}^i \rangle
     - 2 \langle \boldsymbol{v}^i|\boldsymbol{u}^i \rangle + \\
     + 2c\sum_{j>i}\left(\langle \boldsymbol{u}^i|\boldsymbol{u}^j \rangle + \langle \boldsymbol{v}^i|\boldsymbol{v}^j \rangle
       -\langle \boldsymbol{v}^i|\boldsymbol{u}^j \rangle - \langle \boldsymbol{v}^j|\boldsymbol{u}^i \rangle\right)\Bigg]\end{gathered}

Van Rossum-like metrics
-----------------------

In Van Rossum-like metrics, the feature function and the single-unit
kernel are, for :math:`\tau\neq 0`,

.. math::

   \begin{gathered}
   \phi^{\textrm{VR}}_{\tau}(t) = \sqrt{\frac{2}{\tau}}\cdot e^{-t/\tau}\theta(t) \\
   \mathcal{K}^{\textrm{VR}}_{\tau}(t_1,t_2) = \begin{cases}
     1 & \textrm{if } t_1=t_2\\
     e^{-\left\vert t_1-t_2 \right\vert/\tau} & \textrm{if } t_1\neq t_2
   \end{cases}\end{gathered}

where :math:`\theta` is the Heaviside step function (with
:math:`\theta(0)=1`), and we have chosen to normalise
:math:`\phi^{\textrm{VR}}_{\tau}` so that

.. math:: \left\Vert \phi^{\textrm{VR}}_{\tau} \right\Vert_2 = \sqrt{\int\textrm{d }\!t \left[\phi^{\textrm{VR}}_{\tau}(t)\right]^2} = 1 \quad.

In the :math:`\tau\rightarrow 0` limit,

.. math::

   \begin{gathered}
   \phi^{\textrm{VR}}_{0}(t) = \delta(t)\\
   \mathcal{K}^{\textrm{VR}}_{0}(t_1,t_2) = \begin{cases}
     1 & \textrm{if } t_1=t_2\\
     0 & \textrm{if } t_1\neq t_2
   \end{cases}\end{gathered}

In particular, the single-unit inner product now becomes

.. math::

   \langle \boldsymbol{u}|\boldsymbol{v}\rangle = \sum_{n=1}^N\sum_{m=1}^M
     \mathcal{K}^{\textrm{VR}}(u_n,v_m) = \sum_{n=1}^N\sum_{m=1}^M
     e^{-\left\vert u_n-v_m \right\vert/\tau}

Markage formulas
~~~~~~~~~~~~~~~~

For a spike train :math:`\boldsymbol{u}` of length :math:`N` and a time
:math:`t` we define the index
:math:`\tilde{N}\left( \boldsymbol{u}, t\right)`

.. math:: \tilde{N}\left( \boldsymbol{u}, t\right) \doteq \max\{n | u_n < t\}

which we can use to re-write
:math:`\langle \boldsymbol{u}|\boldsymbol{u}\rangle` without the
absolute values:

.. math::

   \begin{split}
       \langle \boldsymbol{u}|\boldsymbol{u}\rangle&= \sum_{n=1}^N \left(
         \sum_{m|v_m<u_n}e^{-(u_n-v_m)/\tau} +
         \sum_{m|v_m>u_n}e^{-(v_m-u_n)/\tau} +
         \sum_{m=1}^M\delta\left(u_n,v_m\right)
       \right)\\
       &= \sum_{n=1}^N \left( \sum_{m|v_m<u_n}e^{-(u_n-v_m)/\tau} +
         \sum_{m|u_m<v_n}e^{-(v_n-u_m)/\tau} +
         \sum_{m=1}^M\delta\left(u_n,v_m\right) \right)\\
       &= \sum_{n=1}^N \left[
         \sum_{m=1}^{\tilde{N}\left( \boldsymbol{v}, u_n\right)}e^{-(u_n-v_m)/\tau} +
         \sum_{m=1}^{\tilde{N}\left( \boldsymbol{u}, v_n\right)}e^{-(v_n-u_m)/\tau} +
         \delta\left(u_n,v_{\tilde{N}\left( \boldsymbol{v}, u_n\right)+1}\right)
       \right]\\
       &= \sum_{n=1}^N \Bigg[ e^{-(u_n-v_{\tilde{N}\left( \boldsymbol{v}, u_n\right)})/\tau}
       \sum_{m=1}^{\tilde{N}\left( \boldsymbol{v}, u_n\right)}e^{-(v_{\tilde{N}\left( \boldsymbol{v}, u_n\right)}-v_m)/\tau} +
       \\
       &\phantom{ = \sum_{n=1}^N } + e^{-(v_n-u_{\tilde{N}\left( \boldsymbol{u}, v_n\right)})/\tau}
       \sum_{m=1}^{\tilde{N}\left( \boldsymbol{u}, v_n\right)}e^{-(u_{\tilde{N}\left( \boldsymbol{u}, v_n\right)}-u_m)/\tau}
       + \\
       &\phantom{ = \sum_{n=1}^N } + \delta\left(u_n,v_{\tilde{N}\left( \boldsymbol{v}, u_n\right)+1}\right) \Bigg]\\
     \end{split}

For a spike train :math:`\boldsymbol{u}` of length :math:`N`, we also
define the the *markage vector* :math:`\boldsymbol{m}`, with the same
length as :math:`\boldsymbol{u}`, through the following recursive
assignment:

.. math::

   \begin{aligned}
     m_1(\boldsymbol{u}) &\doteq 0 \\
     m_n(\boldsymbol{u}) &\doteq \left(m_{n-1} + 1\right) e^{-(u_n - u_{n-1})/\tau}
     \quad \forall n \in \{2,\ldots,N\}\label{eq:markage_definition}\end{aligned}

It is easy to see that

.. math::

   \begin{split}
       m_n(\boldsymbol{u}) &= \sum_{k=1}^{n-1}e^{-(u_n - u_k)/\tau} =
       \left(\sum_{k=1}^{n}e^{-(u_n - u_k)/\tau}\right) - e^{-(u_n - u_n)/\tau}\\
       &= \sum_{k=1}^{n}e^{-(u_n - u_k)/\tau} - 1
     \end{split}

and in particular

.. math::

   \label{eq:markage_sum}
     \sum_{n=1}^{\tilde{N}\left( \boldsymbol{u}, t\right)}e^{-(u_{\tilde{N}\left( \boldsymbol{u}, t\right)}-u_n)/\tau} = 1 + m_{\tilde{N}\left( \boldsymbol{u}, t\right)}(\boldsymbol{u})

With this definition, we get

.. math::

   \label{eq:singleunit_intprod_markage}
     \begin{split}
       \langle \boldsymbol{u}|\boldsymbol{v}\rangle &= \sum_{n=1}^N \Bigg[
       e^{-(u_n-v_{\tilde{N}\left( \boldsymbol{v}, u_n\right)})/\tau} \left(1 + m_{\tilde{N}\left( \boldsymbol{v}, u_n\right)}(\boldsymbol{v})\right)
       + \\
       &\phantom{= \sum_{n=1}^N} + e^{-(v_n-u_{\tilde{N}\left( \boldsymbol{u}, v_n\right)})/\tau}
       \left(1 + m_{\tilde{N}\left( \boldsymbol{u}, v_n\right)}(\boldsymbol{u})\right) + \\
       &\phantom{= \sum_{n=1}^N} + 
       \delta\left(u_n,v_{\tilde{N}\left( \boldsymbol{v}, u_n\right)+1}\right) \Bigg]
     \end{split}

Finally, note that because of the definition of the markage vector

.. math::

   e^{-(u_n-u_{\tilde{N}\left( \boldsymbol{u}, u_n\right)})/\tau} \left(1 +
       m_{\tilde{N}\left( \boldsymbol{u}, u_n\right)}(\boldsymbol{u})\right) = e^{-(u_n-u_{n-1})/\tau}\left(1+m(\boldsymbol{u})\right) = m_{n}(\boldsymbol{u})

so that in particular

.. math::

   \label{eq:singleunit_squarenorm_markage}
     \begin{split}
       \langle \boldsymbol{u}|\boldsymbol{u}\rangle &= \sum_{n=1}^N \left(1 + 2m_{n}(\boldsymbol{u})\right)
     \end{split}

A formula for the efficient computation of multiunit Van Rossum spike
train metrics, used by ``pymuvr``, can then be obtained by opportunely
substituting these expressions for the single-unit scalar products in
the definition of the multiunit distance.
